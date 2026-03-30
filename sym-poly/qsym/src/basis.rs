use std::fmt;

/// Bases for the ring of quasisymmetric functions.
///
/// Bases for the ring of quasisymmetric functions.
///
/// - `Monomial`      -- M_α (quasisymmetric monomial)
/// - `Fundamental`   -- F_α (Gessel fundamental)
/// - `PowerSumPsi`   -- Ψ_α (type 1 power sum)
/// - `PowerSumPhi`   -- Φ_α (type 2 power sum)
///
/// The power sum bases Ψ and Φ are defined in:
/// Ballantine--Daugherty--Hicks--Mason--Niese,
/// *Quasisymmetric Power Sums*, JCTA 2020.
/// <https://doi.org/10.1016/j.jcta.2020.105273>
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum QSymBasis {
    Monomial,
    Fundamental,
    PowerSumPsi,
    PowerSumPhi,
}

impl QSymBasis {
    /// Short symbol used in display.
    pub fn symbol(&self) -> &'static str {
        match self {
            QSymBasis::Monomial => "M",
            QSymBasis::Fundamental => "F",
            QSymBasis::PowerSumPsi => "Ψ",
            QSymBasis::PowerSumPhi => "Φ",
        }
    }
}

impl fmt::Display for QSymBasis {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.symbol())
    }
}
