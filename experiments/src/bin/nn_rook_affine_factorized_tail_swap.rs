//! Probe the path-word involution needed by the affine LGV route.
//!
//! The algebraic check in `nn_rook_lgv_affine_transfer_certificate_match`
//! verifies that
//!
//!   A paths + t-marked Q paths + s-marked shifted Q paths
//!
//! have the right coefficient sequence.  This binary asks for more: if two
//! such path words have inverted Toeplitz endpoints, does the usual first
//! intersection tail swap close back into the same factorized path language?
//!
//! The model is intentionally concrete and conservative.  A generated word has
//! one of the transfer components used in the algebraic check, row events with
//! physical scan intervals after stripping, and a terminal `s`/`t` marker plus
//! optional degree shift.  At a shared select interval we splice the prefix of
//! one word to the suffix/terminal marker of the other word and then test
//! whether the crossed words are present among the originally generated words.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::env;

use experiments::nn_rook_utils::partitions;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
enum QKind {
    Upper,
    Lower,
}

impl QKind {
    fn name(self) -> &'static str {
        match self {
            QKind::Upper => "U",
            QKind::Lower => "L",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
enum Param {
    None,
    S,
    T,
}

#[derive(Clone, Copy)]
struct TermOptions {
    include_a: bool,
    include_s: bool,
    include_t: bool,
}

impl TermOptions {
    fn parse(input: &str) -> Self {
        match input {
            "no_s" => Self {
                include_a: true,
                include_s: false,
                include_t: true,
            },
            "no_t" => Self {
                include_a: true,
                include_s: true,
                include_t: false,
            },
            "q_only" => Self {
                include_a: false,
                include_s: true,
                include_t: true,
            },
            "a_only" => Self {
                include_a: true,
                include_s: false,
                include_t: false,
            },
            _ => Self {
                include_a: true,
                include_s: true,
                include_t: true,
            },
        }
    }
}

impl Param {
    fn token(self) -> Option<&'static str> {
        match self {
            Param::None => None,
            Param::S => Some("s"),
            Param::T => Some("t"),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
enum Component {
    ABase { strip: usize },
    ATail { strip: usize },
    UBase { strip: usize },
    UEnd { strip: usize },
    UReservoir { strip: usize },
    LReservoir { strip: usize },
    LWindow { strip: usize, window_width: usize },
}

#[derive(Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
enum RowEvent {
    Skip {
        bound_phys: usize,
    },
    Select {
        lower_phys: usize,
        upper_phys: usize,
    },
    Reservoir {
        bound_phys: usize,
    },
    TerminalSkip {
        bound_phys: usize,
    },
}

impl RowEvent {
    fn is_select(&self) -> bool {
        matches!(self, RowEvent::Select { .. })
    }

    fn weight_row(&self) -> bool {
        matches!(self, RowEvent::Select { .. } | RowEvent::Reservoir { .. })
    }

    fn interval(&self) -> Option<(usize, usize)> {
        match self {
            RowEvent::Select {
                lower_phys,
                upper_phys,
            } => Some((*lower_phys, *upper_phys)),
            _ => None,
        }
    }

    fn target_bound_phys(&self) -> usize {
        match self {
            RowEvent::Skip { bound_phys }
            | RowEvent::Reservoir { bound_phys }
            | RowEvent::TerminalSkip { bound_phys } => *bound_phys,
            RowEvent::Select { lower_phys, .. } => *lower_phys,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
struct Word {
    component: Component,
    param: Param,
    offset: usize,
    events: Vec<RowEvent>,
}

impl Word {
    fn selected_rows(&self) -> usize {
        self.events.iter().filter(|event| event.is_select()).count()
    }

    fn total_degree(&self) -> usize {
        self.selected_rows() + self.offset
    }

    fn weight_signature(&self) -> BTreeMap<String, usize> {
        let mut out = BTreeMap::new();
        for (row, event) in self.events.iter().enumerate() {
            if event.weight_row() {
                *out.entry(format!("w{row}")).or_insert(0) += 1;
            }
        }
        if let Some(token) = self.param.token() {
            *out.entry(token.to_string()).or_insert(0) += 1;
        }
        out
    }
}

#[derive(Clone, Debug)]
struct AbsPath {
    source: usize,
    sink: usize,
    source_degree: usize,
    sink_degree: usize,
    word: Word,
}

impl AbsPath {
    fn key(&self) -> PathKey {
        PathKey {
            source: self.source,
            sink: self.sink,
            word: self.word.clone(),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
struct PathKey {
    source: usize,
    sink: usize,
    word: Word,
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
struct LoosePathKey {
    source: usize,
    sink: usize,
    param: Param,
    offset: usize,
    events: Vec<RowEvent>,
}

impl AbsPath {
    fn loose_key(&self) -> LoosePathKey {
        LoosePathKey {
            source: self.source,
            sink: self.sink,
            param: self.word.param,
            offset: self.word.offset,
            events: self.word.events.clone(),
        }
    }
}

impl LoosePathKey {
    fn from_word(source: usize, sink: usize, word: &Word) -> Self {
        Self {
            source,
            sink,
            param: word.param,
            offset: word.offset,
            events: word.events.clone(),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
enum Resource {
    Touch {
        row: usize,
        degree: usize,
        phys: usize,
    },
    TargetSkip {
        row: usize,
        degree: usize,
        bound_phys: usize,
    },
    Reservoir {
        row: usize,
        degree: usize,
        bound_phys: usize,
    },
}

#[derive(Clone, Debug)]
enum Candidate {
    Select(SelectSwapSite),
    SuffixAfterRow { row: usize },
}

#[derive(Clone, Debug)]
struct SelectSwapSite {
    row: usize,
    p_interval: (usize, usize),
    q_interval: (usize, usize),
}

#[derive(Default)]
struct Summary {
    tested_pairs: usize,
    no_shared_resource: usize,
    no_swap_candidate: usize,
    product_fail: usize,
    closure_fail: usize,
    closed: usize,
    strict_closed: usize,
    closed_by_select: usize,
    closed_by_suffix: usize,
    first_failure: Option<String>,
    first_no_shared: Option<String>,
    first_no_swap: Option<String>,
    first_closure_failure: Option<String>,
}

fn physical_bound(strip: usize, bound: usize) -> usize {
    strip + bound
}

fn physical_lower(strip: usize, q: usize) -> usize {
    strip + q
}

fn max_stripped_width(eta: &[usize], strip: usize) -> usize {
    eta.iter()
        .map(|&width| width.saturating_sub(strip))
        .max()
        .unwrap_or(0)
}

fn gen_strip_words(
    eta: &[usize],
    h: usize,
    strip: usize,
    component: Component,
    param: Param,
    offset: usize,
) -> Vec<Word> {
    let start_bound = max_stripped_width(eta, strip) + 1;
    let mut out = Vec::new();
    let mut events = Vec::with_capacity(eta.len());

    fn rec(
        eta: &[usize],
        h: usize,
        strip: usize,
        row: usize,
        bound: usize,
        component: &Component,
        param: Param,
        offset: usize,
        events: &mut Vec<RowEvent>,
        out: &mut Vec<Word>,
    ) {
        if row == eta.len() {
            out.push(Word {
                component: component.clone(),
                param,
                offset,
                events: events.clone(),
            });
            return;
        }

        events.push(RowEvent::Skip {
            bound_phys: physical_bound(strip, bound),
        });
        rec(
            eta,
            h,
            strip,
            row + 1,
            bound,
            component,
            param,
            offset,
            events,
            out,
        );
        events.pop();

        let width = eta[row].saturating_sub(strip);
        if row <= h && width > 0 {
            for q in 1..=width.min(bound.saturating_sub(1)) {
                events.push(RowEvent::Select {
                    lower_phys: physical_lower(strip, q),
                    upper_phys: physical_bound(strip, bound),
                });
                rec(
                    eta,
                    h,
                    strip,
                    row + 1,
                    q,
                    component,
                    param,
                    offset,
                    events,
                    out,
                );
                events.pop();
            }
        }
    }

    rec(
        eta,
        h,
        strip,
        0,
        start_bound,
        &component,
        param,
        offset,
        &mut events,
        &mut out,
    );
    out
}

fn gen_reservoir_words(
    eta: &[usize],
    h: usize,
    c: usize,
    component: Component,
    param: Param,
    offset: usize,
) -> Vec<Word> {
    let strip = c + 1;
    let start_bound = max_stripped_width(eta, strip) + 1;
    let mut out = Vec::new();
    let mut events = Vec::with_capacity(eta.len());

    fn rec(
        eta: &[usize],
        h: usize,
        c: usize,
        strip: usize,
        row: usize,
        bound: usize,
        terminal: bool,
        component: &Component,
        param: Param,
        offset: usize,
        events: &mut Vec<RowEvent>,
        out: &mut Vec<Word>,
    ) {
        if row == eta.len() {
            if terminal {
                out.push(Word {
                    component: component.clone(),
                    param,
                    offset,
                    events: events.clone(),
                });
            }
            return;
        }

        if terminal {
            events.push(RowEvent::TerminalSkip {
                bound_phys: physical_bound(strip, bound),
            });
            rec(
                eta,
                h,
                c,
                strip,
                row + 1,
                bound,
                true,
                component,
                param,
                offset,
                events,
                out,
            );
            events.pop();
            return;
        }

        events.push(RowEvent::Skip {
            bound_phys: physical_bound(strip, bound),
        });
        rec(
            eta,
            h,
            c,
            strip,
            row + 1,
            bound,
            false,
            component,
            param,
            offset,
            events,
            out,
        );
        events.pop();

        if eta[row] > c {
            events.push(RowEvent::Reservoir {
                bound_phys: physical_bound(strip, bound),
            });
            rec(
                eta,
                h,
                c,
                strip,
                row + 1,
                bound,
                true,
                component,
                param,
                offset,
                events,
                out,
            );
            events.pop();
        }

        let width = eta[row].saturating_sub(strip);
        if row <= h && width > 0 {
            for q in 1..=width.min(bound.saturating_sub(1)) {
                events.push(RowEvent::Select {
                    lower_phys: physical_lower(strip, q),
                    upper_phys: physical_bound(strip, bound),
                });
                rec(
                    eta,
                    h,
                    c,
                    strip,
                    row + 1,
                    q,
                    false,
                    component,
                    param,
                    offset,
                    events,
                    out,
                );
                events.pop();
            }
        }
    }

    rec(
        eta,
        h,
        c,
        strip,
        0,
        start_bound,
        false,
        &component,
        param,
        offset,
        &mut events,
        &mut out,
    );
    out
}

fn gen_window_words(
    eta: &[usize],
    h: usize,
    c: usize,
    p: usize,
    component: Component,
    param: Param,
    offset: usize,
) -> Vec<Word> {
    let strip = c + 1;
    let window_width = p - c;
    let start_bound = max_stripped_width(eta, strip) + 1;
    let mut out = Vec::new();
    let mut events = Vec::with_capacity(eta.len());

    fn rec(
        eta: &[usize],
        h: usize,
        strip: usize,
        window_width: usize,
        row: usize,
        bound: usize,
        seen: bool,
        component: &Component,
        param: Param,
        offset: usize,
        events: &mut Vec<RowEvent>,
        out: &mut Vec<Word>,
    ) {
        if row == eta.len() {
            if seen {
                out.push(Word {
                    component: component.clone(),
                    param,
                    offset,
                    events: events.clone(),
                });
            }
            return;
        }

        events.push(RowEvent::Skip {
            bound_phys: physical_bound(strip, bound),
        });
        rec(
            eta,
            h,
            strip,
            window_width,
            row + 1,
            bound,
            seen,
            component,
            param,
            offset,
            events,
            out,
        );
        events.pop();

        let width = eta[row].saturating_sub(strip);
        if row <= h && width > 0 {
            for q in 1..=width.min(bound.saturating_sub(1)) {
                events.push(RowEvent::Select {
                    lower_phys: physical_lower(strip, q),
                    upper_phys: physical_bound(strip, bound),
                });
                rec(
                    eta,
                    h,
                    strip,
                    window_width,
                    row + 1,
                    q,
                    seen || q <= window_width,
                    component,
                    param,
                    offset,
                    events,
                    out,
                );
                events.pop();
            }
        }
    }

    rec(
        eta,
        h,
        strip,
        window_width,
        0,
        start_bound,
        false,
        &component,
        param,
        offset,
        &mut events,
        &mut out,
    );
    out
}

fn q_components(eta: &[usize], c: usize, p: usize, h: usize, kind: QKind) -> Vec<Word> {
    match kind {
        QKind::Upper => {
            let mut out = gen_strip_words(
                eta,
                h,
                c + 1,
                Component::UBase { strip: c + 1 },
                Param::None,
                0,
            );
            out.extend(gen_strip_words(
                eta,
                h,
                p + 1,
                Component::UEnd { strip: p + 1 },
                Param::None,
                0,
            ));
            out.extend(gen_reservoir_words(
                eta,
                h,
                c,
                Component::UReservoir { strip: c + 1 },
                Param::None,
                0,
            ));
            out
        }
        QKind::Lower => {
            let mut out = gen_reservoir_words(
                eta,
                h,
                c,
                Component::LReservoir { strip: c + 1 },
                Param::None,
                0,
            );
            out.extend(gen_window_words(
                eta,
                h,
                c,
                p,
                Component::LWindow {
                    strip: c + 1,
                    window_width: p - c,
                },
                Param::None,
                0,
            ));
            out
        }
    }
}

fn factorized_words(
    eta: &[usize],
    c: usize,
    p: usize,
    h: usize,
    kind: QKind,
    terms: TermOptions,
) -> Vec<Word> {
    let mut out = Vec::new();
    if terms.include_a {
        out.extend(gen_strip_words(
            eta,
            h,
            c + 1,
            Component::ABase { strip: c + 1 },
            Param::None,
            0,
        ));
        for strip in (c + 2)..=(p + 1) {
            out.extend(gen_strip_words(
                eta,
                h,
                strip,
                Component::ATail { strip },
                Param::None,
                1,
            ));
        }
    }

    for mut word in q_components(eta, c, p, h, kind) {
        if terms.include_t {
            word.param = Param::T;
            out.push(word.clone());
        }
        if terms.include_s {
            word.param = Param::S;
            word.offset += 1;
            out.push(word);
        }
    }
    out
}

fn all_abs_paths(words: &[Word], rows: &[usize], cols: &[usize]) -> Vec<AbsPath> {
    let mut out = Vec::new();
    for (source, &source_degree) in cols.iter().enumerate() {
        for (sink, &sink_degree) in rows.iter().enumerate() {
            if sink_degree < source_degree {
                continue;
            }
            let diff = sink_degree - source_degree;
            for word in words {
                if word.total_degree() == diff {
                    out.push(AbsPath {
                        source,
                        sink,
                        source_degree,
                        sink_degree,
                        word: word.clone(),
                    });
                }
            }
        }
    }
    out
}

fn resources(path: &AbsPath, global_max_bound: usize) -> BTreeSet<Resource> {
    let mut out = BTreeSet::new();
    let mut selected_before = 0usize;
    for (row, event) in path.word.events.iter().enumerate() {
        let degree = path.source_degree + selected_before;
        match event {
            RowEvent::Skip { bound_phys } | RowEvent::TerminalSkip { bound_phys } => {
                out.insert(Resource::TargetSkip {
                    row,
                    degree,
                    bound_phys: *bound_phys,
                });
            }
            RowEvent::Select {
                lower_phys,
                upper_phys,
            } => {
                for phys in *lower_phys..*upper_phys {
                    out.insert(Resource::Touch { row, degree, phys });
                }
                for bound_phys in 1..=global_max_bound {
                    out.insert(Resource::TargetSkip {
                        row,
                        degree: degree + 1,
                        bound_phys,
                    });
                }
                selected_before += 1;
            }
            RowEvent::Reservoir { bound_phys } => {
                out.insert(Resource::Reservoir {
                    row,
                    degree,
                    bound_phys: *bound_phys,
                });
                out.insert(Resource::TargetSkip {
                    row,
                    degree,
                    bound_phys: *bound_phys,
                });
            }
        }
    }
    out
}

fn first_select_swap_site(p: &AbsPath, q: &AbsPath) -> Option<SelectSwapSite> {
    let mut p_selected = 0usize;
    let mut q_selected = 0usize;
    for row in 0..p.word.events.len() {
        let p_degree = p.source_degree + p_selected;
        let q_degree = q.source_degree + q_selected;
        let p_interval = p.word.events[row].interval();
        let q_interval = q.word.events[row].interval();
        if let (Some((p_lower, p_upper)), Some((q_lower, q_upper))) = (p_interval, q_interval) {
            if p_degree == q_degree && p_lower.max(q_lower) < p_upper.min(q_upper) {
                return Some(SelectSwapSite {
                    row,
                    p_interval: (p_lower, p_upper),
                    q_interval: (q_lower, q_upper),
                });
            }
        }
        if p.word.events[row].is_select() {
            p_selected += 1;
        }
        if q.word.events[row].is_select() {
            q_selected += 1;
        }
    }
    None
}

fn degree_after_row(path: &AbsPath, row: usize) -> usize {
    path.source_degree
        + path
            .word
            .events
            .iter()
            .take(row + 1)
            .filter(|event| event.is_select())
            .count()
}

fn resource_row(resource: &Resource) -> usize {
    match resource {
        Resource::Touch { row, .. }
        | Resource::TargetSkip { row, .. }
        | Resource::Reservoir { row, .. } => *row,
    }
}

fn row_has_shared_resource(
    row: usize,
    left_resources: &BTreeSet<Resource>,
    right_resources: &BTreeSet<Resource>,
) -> bool {
    left_resources
        .intersection(right_resources)
        .any(|resource| resource_row(resource) == row)
}

fn swap_candidates(
    p: &AbsPath,
    q: &AbsPath,
    left_resources: &BTreeSet<Resource>,
    right_resources: &BTreeSet<Resource>,
) -> Vec<Candidate> {
    let mut out = Vec::new();
    if let Some(site) = first_select_swap_site(p, q) {
        out.push(Candidate::Select(site));
    }
    for row in 0..p.word.events.len() {
        let last_row = row + 1 == p.word.events.len();
        let bounds_match =
            p.word.events[row].target_bound_phys() == q.word.events[row].target_bound_phys();
        if degree_after_row(p, row) == degree_after_row(q, row)
            && (last_row || bounds_match)
            && row_has_shared_resource(row, left_resources, right_resources)
        {
            out.push(Candidate::SuffixAfterRow { row });
        }
    }
    out
}

fn combine_signatures(
    a: &BTreeMap<String, usize>,
    b: &BTreeMap<String, usize>,
) -> BTreeMap<String, usize> {
    let mut out = a.clone();
    for (key, value) in b {
        *out.entry(key.clone()).or_insert(0) += value;
    }
    out
}

fn crossed_words(p: &AbsPath, q: &AbsPath, candidate: &Candidate) -> Option<(Word, Word)> {
    match candidate {
        Candidate::Select(site) => {
            let (p_lower, p_upper) = site.p_interval;
            let (q_lower, q_upper) = site.q_interval;
            if q_lower >= p_upper || p_lower >= q_upper {
                return None;
            }

            let mut left_events = Vec::new();
            left_events.extend_from_slice(&p.word.events[..site.row]);
            left_events.push(RowEvent::Select {
                lower_phys: q_lower,
                upper_phys: p_upper,
            });
            left_events.extend_from_slice(&q.word.events[(site.row + 1)..]);

            let mut right_events = Vec::new();
            right_events.extend_from_slice(&q.word.events[..site.row]);
            right_events.push(RowEvent::Select {
                lower_phys: p_lower,
                upper_phys: q_upper,
            });
            right_events.extend_from_slice(&p.word.events[(site.row + 1)..]);

            Some((
                Word {
                    component: q.word.component.clone(),
                    param: q.word.param,
                    offset: q.word.offset,
                    events: left_events,
                },
                Word {
                    component: p.word.component.clone(),
                    param: p.word.param,
                    offset: p.word.offset,
                    events: right_events,
                },
            ))
        }
        Candidate::SuffixAfterRow { row } => {
            let split = row + 1;
            if split == p.word.events.len() {
                return Some((
                    Word {
                        component: p.word.component.clone(),
                        param: q.word.param,
                        offset: q.word.offset,
                        events: p.word.events.clone(),
                    },
                    Word {
                        component: q.word.component.clone(),
                        param: p.word.param,
                        offset: p.word.offset,
                        events: q.word.events.clone(),
                    },
                ));
            }

            let mut left_events = Vec::new();
            left_events.extend_from_slice(&p.word.events[..split]);
            left_events.extend_from_slice(&q.word.events[split..]);

            let mut right_events = Vec::new();
            right_events.extend_from_slice(&q.word.events[..split]);
            right_events.extend_from_slice(&p.word.events[split..]);

            Some((
                Word {
                    component: q.word.component.clone(),
                    param: q.word.param,
                    offset: q.word.offset,
                    events: left_events,
                },
                Word {
                    component: p.word.component.clone(),
                    param: p.word.param,
                    offset: p.word.offset,
                    events: right_events,
                },
            ))
        }
    }
}

fn describe_path(path: &AbsPath) -> String {
    format!(
        "src {}({}) -> sink {}({}), total_degree={}, component={:?}, param={:?}, offset={}, events={:?}",
        path.source,
        path.source_degree,
        path.sink,
        path.sink_degree,
        path.word.total_degree(),
        path.word.component,
        path.word.param,
        path.word.offset,
        path.word.events
    )
}

fn check_packet(
    eta: &[usize],
    c: usize,
    p: usize,
    h: usize,
    kind: QKind,
    terms: TermOptions,
    minor_size: usize,
    summary: &mut Summary,
) {
    let rows: Vec<_> = (0..minor_size).collect();
    let cols = rows.clone();
    let words = factorized_words(eta, c, p, h, kind, terms);
    let paths = all_abs_paths(&words, &rows, &cols);
    let membership: HashSet<PathKey> = paths.iter().map(AbsPath::key).collect();
    let loose_membership: HashSet<LoosePathKey> = paths.iter().map(AbsPath::loose_key).collect();
    let by_entry: HashMap<(usize, usize), Vec<&AbsPath>> = {
        let mut map: HashMap<(usize, usize), Vec<&AbsPath>> = HashMap::new();
        for path in &paths {
            map.entry((path.source, path.sink)).or_default().push(path);
        }
        map
    };
    let global_max_bound = eta.iter().copied().max().unwrap_or(0) + 1;

    for source_left in 0..minor_size {
        for source_right in (source_left + 1)..minor_size {
            for sink_low in 0..minor_size {
                for sink_high in (sink_low + 1)..minor_size {
                    let Some(lefts) = by_entry.get(&(source_left, sink_high)) else {
                        continue;
                    };
                    let Some(rights) = by_entry.get(&(source_right, sink_low)) else {
                        continue;
                    };
                    for &left in lefts {
                        for &right in rights {
                            summary.tested_pairs += 1;
                            let left_resources = resources(left, global_max_bound);
                            let right_resources = resources(right, global_max_bound);
                            if left_resources.is_disjoint(&right_resources) {
                                summary.no_shared_resource += 1;
                                let message = format!(
                                        "no shared resource\n  eta={eta:?} c={c} p={p} h={h} Q={}\n  left={}\n  right={}",
                                        kind.name(),
                                        describe_path(left),
                                        describe_path(right)
                                    );
                                if summary.first_no_shared.is_none() {
                                    summary.first_no_shared = Some(message.clone());
                                }
                                if summary.first_failure.is_none() {
                                    summary.first_failure = Some(message);
                                }
                                continue;
                            }

                            let candidates =
                                swap_candidates(left, right, &left_resources, &right_resources);
                            if candidates.is_empty() {
                                summary.no_swap_candidate += 1;
                                if summary.first_failure.is_none() {
                                    let shared = left_resources
                                        .intersection(&right_resources)
                                        .next()
                                        .cloned();
                                    let message = format!(
                                        "shared resource but no modeled swap candidate\n  eta={eta:?} c={c} p={p} h={h} Q={}\n  shared={shared:?}\n  left={}\n  right={}",
                                        kind.name(),
                                        describe_path(left),
                                        describe_path(right)
                                    );
                                    if summary.first_no_swap.is_none() {
                                        summary.first_no_swap = Some(message.clone());
                                    }
                                    summary.first_failure = Some(message);
                                }
                                continue;
                            };

                            let original_sig = combine_signatures(
                                &left.word.weight_signature(),
                                &right.word.weight_signature(),
                            );
                            let mut saw_product_failure = false;
                            let mut first_closed_candidate = None;
                            let mut first_closure_failure = None;
                            for candidate in &candidates {
                                let Some((cross_left, cross_right)) =
                                    crossed_words(left, right, candidate)
                                else {
                                    continue;
                                };
                                let crossed_sig = combine_signatures(
                                    &cross_left.weight_signature(),
                                    &cross_right.weight_signature(),
                                );
                                if original_sig != crossed_sig {
                                    saw_product_failure = true;
                                    if summary.first_failure.is_none() {
                                        summary.first_failure = Some(format!(
                                            "product signature changed\n  eta={eta:?} c={c} p={p} h={h} Q={}\n  candidate={candidate:?}\n  left={}\n  right={}\n  crossed_left={cross_left:?}\n  crossed_right={cross_right:?}\n  original_sig={original_sig:?}\n  crossed_sig={crossed_sig:?}",
                                            kind.name(),
                                            describe_path(left),
                                            describe_path(right)
                                        ));
                                    }
                                    continue;
                                }

                                let cross_left_key = PathKey {
                                    source: source_left,
                                    sink: sink_low,
                                    word: cross_left.clone(),
                                };
                                let cross_right_key = PathKey {
                                    source: source_right,
                                    sink: sink_high,
                                    word: cross_right.clone(),
                                };
                                let left_present = membership.contains(&cross_left_key);
                                let right_present = membership.contains(&cross_right_key);
                                let cross_left_loose =
                                    LoosePathKey::from_word(source_left, sink_low, &cross_left);
                                let cross_right_loose =
                                    LoosePathKey::from_word(source_right, sink_high, &cross_right);
                                let left_loose_present =
                                    loose_membership.contains(&cross_left_loose);
                                let right_loose_present =
                                    loose_membership.contains(&cross_right_loose);
                                if left_loose_present && right_loose_present {
                                    first_closed_candidate =
                                        Some((candidate.clone(), left_present && right_present));
                                    break;
                                }
                                if first_closure_failure.is_none() {
                                    first_closure_failure = Some((
                                        candidate.clone(),
                                        cross_left,
                                        cross_right,
                                        left_loose_present,
                                        right_loose_present,
                                    ));
                                }
                            }

                            if let Some((candidate, strictly_present)) = first_closed_candidate {
                                summary.closed += 1;
                                if strictly_present {
                                    summary.strict_closed += 1;
                                }
                                match candidate {
                                    Candidate::Select(_) => summary.closed_by_select += 1,
                                    Candidate::SuffixAfterRow { .. } => {
                                        summary.closed_by_suffix += 1
                                    }
                                }
                            } else if saw_product_failure && first_closure_failure.is_none() {
                                summary.product_fail += 1;
                            } else {
                                summary.closure_fail += 1;
                                if let Some((
                                    candidate,
                                    cross_left,
                                    cross_right,
                                    left_present,
                                    right_present,
                                )) = first_closure_failure
                                {
                                    let message = format!(
                                            "modeled swaps are weight-preserving but not closed in generated path language\n  eta={eta:?} c={c} p={p} h={h} Q={q_name}\n  first_candidate={candidate:?}\n  left={left_desc}\n  right={right_desc}\n  crossed_left source={source_left} sink={sink_low} word={cross_left:?} loose_present={left_present}\n  crossed_right source={source_right} sink={sink_high} word={cross_right:?} loose_present={right_present}",
                                            q_name = kind.name(),
                                            left_desc = describe_path(left),
                                            right_desc = describe_path(right),
                                        );
                                    if summary.first_closure_failure.is_none() {
                                        summary.first_closure_failure = Some(message.clone());
                                    }
                                    if summary.first_failure.is_none() {
                                        summary.first_failure = Some(message);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

fn boundary_packets(eta: &[usize]) -> Vec<(usize, usize, usize)> {
    let row_add_width = eta[eta.len() - 1];
    let mut out = Vec::new();
    for c in 0..row_add_width {
        for p in (c + 1)..row_add_width {
            for h in 0..eta.len() {
                out.push((c, p, h));
            }
        }
    }
    out
}

fn main() {
    let max_n = env::args()
        .nth(1)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(6);
    let minor_size = env::args()
        .nth(2)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(3);
    let kind_filter = env::args().nth(3).unwrap_or_else(|| "all".to_string());
    let kinds: Vec<QKind> = match kind_filter.as_str() {
        "U" | "u" | "upper" => vec![QKind::Upper],
        "L" | "l" | "lower" => vec![QKind::Lower],
        _ => vec![QKind::Upper, QKind::Lower],
    };
    let term_filter = env::args().nth(4).unwrap_or_else(|| "all".to_string());
    let terms = TermOptions::parse(&term_filter);

    println!("=== Affine factorized path-word tail-swap probe ===");
    println!(
        "Checking |eta| <= {max_n}, minor size {minor_size}, Q filter {kind_filter}, terms {term_filter}.\n"
    );

    let mut packets_checked = 0usize;
    let mut summary = Summary::default();
    for n in 1..=max_n {
        for eta in partitions(n) {
            if eta.is_empty() || eta[eta.len() - 1] < 2 || eta.len() >= 12 {
                continue;
            }
            for (c, p, h) in boundary_packets(&eta) {
                packets_checked += 1;
                for &kind in &kinds {
                    check_packet(&eta, c, p, h, kind, terms, minor_size, &mut summary);
                }
            }
        }
    }

    println!("packets checked: {packets_checked}");
    println!("inverted path pairs tested: {}", summary.tested_pairs);
    println!(
        "closed swaps with component reinterpretation: {}",
        summary.closed
    );
    println!("  strictly same component label: {}", summary.strict_closed);
    println!(
        "  closed by select/select splice: {}",
        summary.closed_by_select
    );
    println!("  closed by suffix swap: {}", summary.closed_by_suffix);
    println!("no shared resource: {}", summary.no_shared_resource);
    println!("shared but no modeled swap: {}", summary.no_swap_candidate);
    println!("product failures: {}", summary.product_fail);
    println!("closure failures: {}", summary.closure_fail);
    if let Some(failure) = summary.first_failure {
        println!("\nfirst issue:\n{failure}");
    } else {
        println!("\nAll tested factorized tail swaps closed.");
    }
    if let Some(failure) = summary.first_no_shared {
        println!("\nfirst no-shared-resource issue:\n{failure}");
    }
    if let Some(failure) = summary.first_no_swap {
        println!("\nfirst no-modeled-swap issue:\n{failure}");
    }
    if let Some(failure) = summary.first_closure_failure {
        println!("\nfirst closure issue:\n{failure}");
    }
}
