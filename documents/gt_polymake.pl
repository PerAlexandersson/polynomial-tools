#!/usr/bin/env polymake --script
#
# Construct the GT polytope GT(lambda, w) and compute its Ehrhart polynomial
# via polymake. Used for benchmarking against our DP+reciprocity method.
#
# Usage: polymake --script gt_polymake.pl
#
# The GT polytope for shape lambda = (l_1, ..., l_n) and weight w = (w_1, ..., w_k)
# is the set of all triangular arrays
#
#   a^(0)_1  a^(0)_2  ...  a^(0)_n      = lambda
#   a^(1)_1  a^(1)_2  ...  a^(1)_n
#   ...
#   a^(k)_1  a^(k)_2  ...  a^(k)_n      = (0, ..., 0)
#
# subject to interlacing: a^(r)_j >= a^(r+1)_j >= a^(r)_{j+1}
# and row sum constraints: sum_j a^(r)_j = sum_{i>r} w_i (cumulative from bottom).
#
# The free variables are the INTERIOR levels a^(1), ..., a^(k-1).

use application "polytope";

# ---- Configuration: from environment or defaults ----
# Usage: GT_LAMBDA=3,2,1 GT_WEIGHT=1,1,1,1,1,1 polymake --script gt_polymake.pl
my @lambda = defined $ENV{GT_LAMBDA} ? split(/,/, $ENV{GT_LAMBDA}) : (3, 2, 1);
my @weight = defined $ENV{GT_WEIGHT} ? split(/,/, $ENV{GT_WEIGHT}) : (1, 1, 1, 1, 1, 1);
# -----------------------------------

my $n = scalar @lambda;   # number of parts
my $k = scalar @weight;   # number of levels (= length of weight)

# Cumulative weight sums from the top: row_sum[r] = sum of w[0..r]
# Level 0 has row sum = |lambda|, level k has row sum = 0.
# Level r has row sum = |lambda| - sum(w[0..r-1]).
my @cumsum;
my $total = 0;
for my $i (0 .. $k-1) { $total += $weight[$i]; push @cumsum, $total; }

# Free variables: levels 1, ..., k-1 (level 0 = lambda, level k = (0,...,0))
# Variable index: level r (1-indexed from 1 to k-1), position j (0-indexed from 0 to n-1)
# var(r, j) = (r-1)*n + j   for r = 1..k-1, j = 0..n-1
my $num_vars = ($k - 1) * $n;

if ($num_vars == 0) {
    print "Trivial case: 0 free variables (point polytope)\n";
    print "Ehrhart polynomial: 1\n";
    exit(0);
}

sub var_idx {
    my ($r, $j) = @_;  # r = 1..k-1, j = 0..n-1
    return ($r - 1) * $n + $j;
}

# Level values: level 0 = lambda, level k = zeros
sub level_val {
    my ($r, $j) = @_;
    if ($r == 0) { return $j < $n ? $lambda[$j] : 0; }
    if ($r == $k) { return 0; }
    return undef;  # free variable
}

# Build inequality system: A x >= b  (polymake convention: first column is constant)
# Each inequality is [b, -a_1, -a_2, ...] meaning b + (-a_1)*x_1 + ... >= 0
# i.e., a_1*x_1 + ... <= b  becomes  b - a_1*x_1 - ... >= 0

my @ineqs;
my @eqs;

for my $r (0 .. $k-1) {
    for my $j (0 .. $n-1) {
        # Interlacing: a^(r)_j >= a^(r+1)_j
        #   If both free: x[r,j] - x[r+1,j] >= 0
        #   If r=0 (fixed): lambda[j] - x[1,j] >= 0  => lambda[j] >= x[1,j]
        #   If r+1=k (fixed): x[r,j] >= 0

        my @row = (0) x ($num_vars + 1);
        my $lhs_fixed = level_val($r, $j);
        my $rhs_fixed = level_val($r + 1, $j);

        if (defined $lhs_fixed && defined $rhs_fixed) {
            # Both fixed: lambda[j] >= 0 (trivially true), skip
            next;
        } elsif (defined $lhs_fixed) {
            # lambda[j] >= x[r+1, j]  => lambda[j] - x[r+1,j] >= 0
            $row[0] = $lhs_fixed;
            $row[var_idx($r + 1, $j) + 1] = -1;
        } elsif (defined $rhs_fixed) {
            # x[r, j] >= rhs_fixed  => x[r,j] - rhs_fixed >= 0
            $row[0] = -$rhs_fixed;
            $row[var_idx($r, $j) + 1] = 1;
        } else {
            # x[r,j] >= x[r+1,j]  => x[r,j] - x[r+1,j] >= 0
            $row[var_idx($r, $j) + 1] = 1;
            $row[var_idx($r + 1, $j) + 1] = -1;
        }
        push @ineqs, \@row;

        # Interlacing: a^(r+1)_j >= a^(r)_{j+1}
        #   (only if j+1 < n)
        if ($j + 1 < $n) {
            my @row2 = (0) x ($num_vars + 1);
            my $rhs2_fixed = level_val($r + 1, $j);
            my $lhs2_fixed = level_val($r, $j + 1);

            if (defined $rhs2_fixed && defined $lhs2_fixed) {
                next;
            } elsif (defined $rhs2_fixed) {
                # rhs2_fixed >= x[r, j+1]
                # Wait, we want a^(r+1)_j >= a^(r)_{j+1}
                # If a^(r+1)_j is fixed: rhs2_fixed >= a^(r)_{j+1}
                # If a^(r)_{j+1} is free: rhs2_fixed - x[r, j+1] >= 0
                if ($r == 0) {
                    # a^(0)_{j+1} = lambda[j+1], a^(1)_j = fixed? No, r+1=1 means a^(1)_j.
                    # Actually rhs2_fixed = level_val(r+1, j). If r+1 = k, it's 0.
                    # And lhs2_fixed = level_val(r, j+1).
                    # We want a^(r+1)_j >= a^(r)_{j+1}
                    # rhs2 = a^(r+1)_j, lhs2 = a^(r)_{j+1}
                    # rhs2_fixed >= lhs2 (which is free) => rhs2_fixed - x >= 0
                    # But wait, lhs2_fixed is defined here...
                }
                # rhs2_fixed >= x[r, j+1] (if r > 0 and r < k)
                # a^(r+1)_j is fixed, a^(r)_{j+1} is free
                $row2[0] = $rhs2_fixed;
                $row2[var_idx($r, $j + 1) + 1] = -1 if $r >= 1 && $r <= $k-1;
            } elsif (defined $lhs2_fixed) {
                # x[r+1, j] >= lhs2_fixed
                $row2[0] = -$lhs2_fixed;
                $row2[var_idx($r + 1, $j) + 1] = 1 if $r + 1 >= 1 && $r + 1 <= $k-1;
            } else {
                # x[r+1, j] >= x[r, j+1]  => x[r+1,j] - x[r,j+1] >= 0
                $row2[var_idx($r + 1, $j) + 1] = 1;
                $row2[var_idx($r, $j + 1) + 1] = -1;
            }
            push @ineqs, \@row2;
        }
    }
}

# Row sum equality constraints: sum_j a^(r)_j = |lambda| - cumsum[r-1]
# For free levels r = 1..k-1:
for my $r (1 .. $k-1) {
    my $target = $total - $cumsum[$r - 1];  # row sum at level r
    my @row = (0) x ($num_vars + 1);
    $row[0] = -$target;  # equality: sum x[r,j] - target = 0
    for my $j (0 .. $n-1) {
        $row[var_idx($r, $j) + 1] = 1;
    }
    push @eqs, \@row;
}

# Construct polytope
my $ineq_matrix = new Matrix<Rational>(\@ineqs);
my $eq_matrix = scalar @eqs > 0 ? new Matrix<Rational>(\@eqs) : undef;

my $p;
if (defined $eq_matrix) {
    $p = new Polytope(INEQUALITIES => $ineq_matrix, EQUATIONS => $eq_matrix);
} else {
    $p = new Polytope(INEQUALITIES => $ineq_matrix);
}

print "lambda = (", join(",", @lambda), ")\n";
print "weight = (", join(",", @weight), ")\n";
print "Dimension: ", $p->DIM, "\n";
print "Vertices: ", $p->N_VERTICES, "\n";

# GT polytopes are rational (not lattice) in general.
# Compute h*-polynomial via Normaliz, then convert to Ehrhart.
use Benchmark qw(:all);
my $t0 = Benchmark->new;
my $hstar = $p->H_STAR_VECTOR;
my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);

print "h*-vector: $hstar\n";
print "Time: ", timestr($td), "\n";
