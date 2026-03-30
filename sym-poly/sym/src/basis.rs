use std::fmt;

/// The six classical bases for the ring of symmetric functions.
///
/// Naming follows SymmetricFunctions.m:
/// - `Monomial`     -- m_λ
/// - `Elementary`   -- e_λ
/// - `CompleteH`    -- h_λ
/// - `PowerSum`     -- p_λ
/// - `Schur`        -- s_λ
/// - `Forgotten`    -- f_λ
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Basis {
    Monomial,
    Elementary,
    CompleteH,
    PowerSum,
    Schur,
    Forgotten,
}

impl Basis {
    /// Short symbol used in display (matches Mathematica notation).
    pub fn symbol(&self) -> &'static str {
        match self {
            Basis::Monomial => "m",
            Basis::Elementary => "e",
            Basis::CompleteH => "h",
            Basis::PowerSum => "p",
            Basis::Schur => "s",
            Basis::Forgotten => "f",
        }
    }

    /// Whether products of basis elements of this type are just partition joins.
    /// True for e, h, p (multiplicative bases).
    pub fn is_multiplicative(&self) -> bool {
        matches!(self, Basis::Elementary | Basis::CompleteH | Basis::PowerSum)
    }

    /// The basis related by the omega involution.
    pub fn omega_dual(&self) -> Basis {
        match self {
            Basis::Monomial => Basis::Forgotten,
            Basis::Forgotten => Basis::Monomial,
            Basis::Elementary => Basis::CompleteH,
            Basis::CompleteH => Basis::Elementary,
            Basis::PowerSum => Basis::PowerSum,
            Basis::Schur => Basis::Schur,
        }
    }
}

impl fmt::Display for Basis {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.symbol())
    }
}
