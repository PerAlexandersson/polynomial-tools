use polynomial_tools::recurrence::{
    format_rational_coeff, parse_rational_coeff, BigRational, RecurrenceJson,
};
use std::fs;
use std::path::{Path, PathBuf};

const BASE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/fixtures/recurrence-benchmarks"
);

#[derive(Debug)]
struct FixtureManifestRow {
    slug: String,
    rows_file: String,
    json_file: String,
}

fn parse_manifest() -> Vec<FixtureManifestRow> {
    let manifest = fs::read_to_string(Path::new(BASE).join("manifest.tsv"))
        .expect("read recurrence benchmark manifest");
    manifest
        .lines()
        .skip(1)
        .map(|line| {
            let cols: Vec<&str> = line.split('\t').collect();
            assert_eq!(cols.len(), 7, "manifest row should have 7 columns: {line}");
            FixtureManifestRow {
                slug: cols[0].to_string(),
                rows_file: cols[5].to_string(),
                json_file: cols[6].to_string(),
            }
        })
        .collect()
}

fn parse_row(line: &str) -> Vec<BigRational> {
    line.split(',')
        .map(|coeff| parse_rational_coeff(coeff.trim()).expect("parse rational coefficient"))
        .collect()
}

fn read_rows(path: PathBuf) -> Vec<Vec<BigRational>> {
    fs::read_to_string(path)
        .expect("read coefficient rows")
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(parse_row)
        .collect()
}

fn format_rows(rows: &[Vec<BigRational>]) -> Vec<String> {
    rows.iter()
        .map(|row| {
            row.iter()
                .map(format_rational_coeff)
                .collect::<Vec<_>>()
                .join(", ")
        })
        .collect()
}

#[test]
fn recurrence_json_fixtures_regenerate_raw_rows() {
    let fixtures = parse_manifest();
    assert_eq!(fixtures.len(), 20, "expected the benchmark fixture suite");

    for fixture in fixtures {
        let expected = read_rows(Path::new(BASE).join(&fixture.rows_file));
        let json_text =
            fs::read_to_string(Path::new(BASE).join(&fixture.json_file)).expect("read JSON file");
        let recurrence_json: RecurrenceJson =
            serde_json::from_str(&json_text).expect("parse recurrence JSON");
        let (recurrence, first_index, initial_polys) = recurrence_json
            .to_recurrence_parts()
            .expect("convert recurrence JSON");
        let generated = recurrence
            .generate_rows_rational(&initial_polys, first_index, expected.len())
            .unwrap_or_else(|err| panic!("generate rows for {}: {err}", fixture.slug));
        assert_eq!(
            format_rows(&generated),
            format_rows(&expected),
            "JSON recurrence should regenerate {}",
            fixture.slug
        );

        let extended = recurrence
            .generate_rows_rational(&initial_polys, first_index, expected.len() + 5)
            .unwrap_or_else(|err| panic!("extend rows for {}: {err}", fixture.slug));
        assert_eq!(extended.len(), expected.len() + 5);
        assert_eq!(
            format_rows(&extended[..expected.len()]),
            format_rows(&expected),
            "extended JSON recurrence should preserve fixture prefix {}",
            fixture.slug
        );
    }
}
