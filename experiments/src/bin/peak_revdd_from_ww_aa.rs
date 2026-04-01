//! Test: does reversed DD at λ+ follow from forward WW + reversed AA at λ+?
//!
//! We have: W_{k-1}^+ ≪ W_k^+ (forward WW, proved)
//!          A_k^+ ≪ A_{k-1}^+ (reversed AA, proved)
//! Need:    S_k ≪ S_{k-1} (reversed DD)
//!
//! Key identities:
//!   A_k^+ = S_k + T_k,  W_k^+ = tS_k + T_k
//!   A_k^+ = W_k^+ - (t-1)S_k
//!   So S_k = (W_k^+ - A_k^+)/(t-1)  (as a polynomial identity)
//!
//! From A_k^+ ≪ A_{k-1}^+ and W_{k-1}^+ ≪ W_k^+:
//!   maybe S_k ≪ S_{k-1} follows from a "sandwich" argument?
//!
//! Also test: A_k^+ ≪ W_k^+ (diagonal AW, proved)
//! and W_k^+ = A_k^+ + (t-1)S_k (shift lemma structure)
//! So if S_k is "small" relative to A_k^+, interlacing is maintained.
//!
//! Alternative: S_k ≪ S_{k-1} iff S_{k-1} ≪ tS_k (Wagner)
//! tS_k = W_k^+ - T_k. And S_{k-1} = A_{k-1}^+ - T_{k-1}.
//! So S_{k-1} ≪ tS_k iff (A_{k-1}^+ - T_{k-1}) ≪ (W_k^+ - T_k).
//! Since T_{k-1} = T_k + W_{k-1}:
//! S_{k-1} = A_{k-1}^+ - T_k - W_{k-1}
//! tS_k = W_k^+ - T_k
//! So S_{k-1} ≪ tS_k becomes:
//! (A_{k-1}^+ - T_k - W_{k-1}) ≪ (W_k^+ - T_k)
//! i.e. (A_{k-1}^+ - W_{k-1}) ≪ W_k^+
//!
//! Now A_{k-1}^+ - W_{k-1} = (S_{k-1} + T_{k-1}) - W_{k-1}
//!   = S_{k-1} + T_{k-1} - W_{k-1} = S_{k-1} + T_k  (since T_{k-1} = T_k + W_{k-1})
//! Wait: T_{k-1} - W_{k-1} = T_k. So A_{k-1}^+ - W_{k-1} = S_{k-1} + T_k = A_k^+!
//!
//! So S_{k-1} ≪ tS_k iff A_k^+ ≪ W_k^+ !!!
//!
//! And A_k^+ ≪ W_k^+ is exactly the DIAGONAL AW at λ+, which IS PROVED (Step 3, k=l)!
//!
//! Therefore: reversed DD at λ+ follows from diagonal AW at λ+!

fn main() {
    println!("=== Key algebraic identity ===");
    println!();
    println!("S_{{k-1}} ≪ tS_k  (Wagner form of reversed DD)");
    println!("  ⟺  (A_{{k-1}}^+ - T_k - W_{{k-1}}) ≪ (W_k^+ - T_k)");
    println!("  ⟺  (S_{{k-1}} + T_{{k-1}} - W_{{k-1}}) ≪ W_k^+");
    println!("  ⟺  (S_{{k-1}} + T_k) ≪ W_k^+        [since T_{{k-1}} - W_{{k-1}} = T_k]");
    println!("  ⟺  A_k^+ ≪ W_k^+                     [since S_{{k-1}} + A_{{k-1}} = S_k, S_k + T_k = A_k^+]");
    println!();
    println!("Wait: S_{{k-1}} + T_k ≠ A_k^+ in general.");
    println!("A_k^+ = S_k + T_k = (S_{{k-1}} + A_{{k-1}}) + T_k");
    println!("So S_{{k-1}} + T_k = A_k^+ - A_{{k-1}}");
    println!();
    println!("Let me redo: S_{{k-1}} ≪ tS_k ⟺ ?");
    println!("tS_k = W_k^+ - T_k");
    println!("S_{{k-1}} = A_{{k-1}}^+ - T_{{k-1}} = A_{{k-1}}^+ - T_k - W_{{k-1}}");
    println!();
    println!("So: (A_{{k-1}}^+ - T_k - W_{{k-1}}) ≪ (W_k^+ - T_k)");
    println!("This is NOT the same as A_k^+ ≪ W_k^+.");
    println!();
    println!("But let's check: what IS A_{{k-1}}^+ - W_{{k-1}}?");
    println!("A_{{k-1}}^+ = S_{{k-1}} + T_{{k-1}} = S_{{k-1}} + T_k + W_{{k-1}}");
    println!("So A_{{k-1}}^+ - W_{{k-1}} = S_{{k-1}} + T_k");
    println!();
    println!("And A_k^+ = S_k + T_k = S_{{k-1}} + A_{{k-1}} + T_k");
    println!("So A_{{k-1}}^+ - W_{{k-1}} = A_k^+ - A_{{k-1}}");
    println!();
    println!("The condition becomes:");
    println!("  (A_k^+ - A_{{k-1}}) ≪ W_k^+");
    println!("  i.e., (A_k^+ - A_{{k-1}}) ≪ A_k^+ + (t-1)S_k");
    println!();
    println!("Hmm, this is (A_k^+ - A_{{k-1}}) ≪ W_k^+.");
    println!("Since A_{{k-1}} has nonneg coeffs, A_k^+ - A_{{k-1}} ≪ A_k^+ trivially?");
    println!("NO — subtracting doesn't preserve interlacing.");
    println!();
    println!("But: A_k^+ - A_{{k-1}} = S_{{k-1}} + T_k.");
    println!("And W_k^+ = tS_k + T_k.");
    println!("So need: (S_{{k-1}} + T_k) ≪ (tS_k + T_k).");
    println!("By cone: S_{{k-1}} ≪ tS_k AND T_k ≪ tS_k AND S_{{k-1}} ≪ T_k AND T_k ≪ T_k.");
    println!("Wait, cone on (tS_k + T_k): need S_{{k-1}} ≪ tS_k, S_{{k-1}} ≪ T_k, T_k ≪ tS_k, T_k ≪ T_k.");
    println!("S_{{k-1}} ≪ tS_k is what we're trying to prove (Wagner of reversed DD)!");
    println!("CIRCULAR.");
    println!();
    println!("DIFFERENT approach: (S_{{k-1}} + T_k) ≪ (tS_k + T_k)");
    println!("= S_{{k-1}} + T_k ≪ t(S_{{k-1}} + A_{{k-1}}) + T_k");
    println!("= S_{{k-1}} + T_k ≪ tS_{{k-1}} + tA_{{k-1}} + T_k");
    println!("= (S_{{k-1}} + T_k) ≪ (S_{{k-1}} + T_k) + (t-1)S_{{k-1}} + tA_{{k-1}}");
    println!();
    println!("Hmm, this is f ≪ f + (t-1)S_{{k-1}} + tA_{{k-1}}.");
    println!("The extra terms are (t-1)S_{{k-1}} + tA_{{k-1}} which has a factor of t...");
    println!("Not a clean (t-1) form.");
    println!();
    println!("YET ANOTHER approach:");
    println!("W_k^+ - (S_{{k-1}} + T_k) = tS_k + T_k - S_{{k-1}} - T_k = tS_k - S_{{k-1}}");
    println!("= t(S_{{k-1}} + A_{{k-1}}) - S_{{k-1}} = (t-1)S_{{k-1}} + tA_{{k-1}}");
    println!();
    println!("So W_k^+ = (S_{{k-1}} + T_k) + (t-1)S_{{k-1}} + tA_{{k-1}}");
    println!("= (S_{{k-1}} + T_k) + (t-1)(S_{{k-1}} + A_{{k-1}}) + A_{{k-1}}");
    println!("= (S_{{k-1}} + T_k) + (t-1)S_k + A_{{k-1}}");
    println!("Hmm.");
    println!();
    println!("One more try: reversed DD for consecutive indices.");
    println!("S_k ≪ S_{{k-1}} ⟺ S_{{k-1}} ≪ tS_k (Wagner).");
    println!();
    println!("tS_k = tS_{{k-1}} + tA_{{k-1}}.");
    println!("Shift lemma on tS_{{k-1}}: tS_k = tS_{{k-1}} + t·A_{{k-1}}.");
    println!("NOT (t-1) form. Can't use shift lemma.");
    println!();
    println!("BUT: what if we write tS_k = S_{{k-1}} + (t-1)S_{{k-1}} + tA_{{k-1}}");
    println!("= S_{{k-1}} + (t-1)(S_{{k-1}} + tA_{{k-1}}/(t-1))... doesn't factor.");
    println!();
    println!("FINAL IDEA: use W_k^+ directly.");
    println!("W_k^+ = tS_k + T_k. And A_k^+ = S_k + T_k.");
    println!("W_k^+ = A_k^+ + (t-1)S_k.");
    println!("Shift lemma: S_k ≪ A_k^+ ⟹ A_k^+ ≪ W_k^+ (diagonal AW). ✓");
    println!();
    println!("Now: A_k^+ ≪ W_k^+ means roots of A_k^+ are more negative than W_k^+.");
    println!("And: A_{{k-1}}^+ ≪ W_{{k-1}}^+ (also diagonal AW). ✓");
    println!("And: W_{{k-1}}^+ ≪ W_k^+ (forward WW). ✓");
    println!("And: A_k^+ ≪ A_{{k-1}}^+ (reversed AA). ✓");
    println!();
    println!("So: A_k^+ ≪ A_{{k-1}}^+ ≪ W_{{k-1}}^+ ≪ W_k^+");
    println!("(a chain, though transitivity doesn't hold in general)");
    println!();
    println!("For reversed DD: S_k ≪ S_{{k-1}}.");
    println!("S_k = A_k^+ - T_k and S_{{k-1}} = A_{{k-1}}^+ - T_{{k-1}}.");
    println!("We know A_k^+ ≪ A_{{k-1}}^+ and T_{{k-1}} = T_k + W_{{k-1}}.");
    println!("So S_{{k-1}} = A_{{k-1}}^+ - T_k - W_{{k-1}} and S_k = A_k^+ - T_k.");
    println!("If we subtract T_k from both sides of A_k^+ ≪ A_{{k-1}}^+...");
    println!("INTERLACING IS NOT PRESERVED UNDER SUBTRACTION.");
}
