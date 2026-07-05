use std::fmt;

/// Bases for the ring of quasisymmetric functions.
///
/// - `Monomial`      -- M_α (quasisymmetric monomial)
/// - `Fundamental`   -- F_α (Gessel fundamental)
/// - `QuasisymmetricSchur` -- S_α (Haglund--Luoto--Mason--van Willigenburg)
/// - `DualImmaculate` -- S*_α (dual immaculate)
/// - `ExtendedSchur` -- E_α (extended Schur)
/// - `RowStrictExtendedSchur` -- row-strict extended Schur
/// - `FlippedExtendedSchur` -- flipped extended Schur
/// - `BackwardExtendedSchur` -- backward extended Schur
/// - `PowerSumPsi`   -- Ψ_α (type 1 power sum)
/// - `PowerSumPhi`   -- Φ_α (type 2 power sum)
/// - `CombinatorialPowerSum` -- p_α (P-partition combinatorial power sum)
/// - `ReverseCombinatorialPowerSum` -- p^r_α (reverse combinatorial power sum)
///
/// The power sum bases Ψ and Φ are defined in:
/// Ballantine--Daugherty--Hicks--Mason--Niese,
/// *Quasisymmetric Power Sums*, JCTA 2020.
/// <https://doi.org/10.1016/j.jcta.2020.105273>
///
/// The combinatorial power sum bases are defined in:
/// Aliniaeifard--Wang--van Willigenburg,
/// *P-partition power sums*, European J. Combin. 2023.
/// <https://doi.org/10.1016/j.ejc.2023.103688>
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum QSymBasis {
    Monomial,
    Fundamental,
    QuasisymmetricSchur,
    DualImmaculate,
    ExtendedSchur,
    RowStrictExtendedSchur,
    FlippedExtendedSchur,
    BackwardExtendedSchur,
    PowerSumPsi,
    PowerSumPhi,
    CombinatorialPowerSum,
    ReverseCombinatorialPowerSum,
}

impl QSymBasis {
    /// Short symbol used in display.
    pub fn symbol(&self) -> &'static str {
        match self {
            QSymBasis::Monomial => "M",
            QSymBasis::Fundamental => "F",
            QSymBasis::QuasisymmetricSchur => "QS",
            QSymBasis::DualImmaculate => "S*",
            QSymBasis::ExtendedSchur => "E",
            QSymBasis::RowStrictExtendedSchur => "RSE",
            QSymBasis::FlippedExtendedSchur => "FE",
            QSymBasis::BackwardExtendedSchur => "BE",
            QSymBasis::PowerSumPsi => "Ψ",
            QSymBasis::PowerSumPhi => "Φ",
            QSymBasis::CombinatorialPowerSum => "p",
            QSymBasis::ReverseCombinatorialPowerSum => "pr",
        }
    }
}

impl fmt::Display for QSymBasis {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.symbol())
    }
}
