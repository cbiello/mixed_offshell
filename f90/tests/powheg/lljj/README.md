# Four-quark channels at O(alpha_s alpha^3): mixed_ll against POWHEG `NSC-lnujjMIX`

l l' + two quarks, tree interference between the O(g_s^2 e^2) and the O(e^4) amplitudes: the hard events of the
double-real sectors `ns_qqb`, `ns_qq` (singular) and `ns_qqb_w`, `ns_qqp_w`, `ns_qqpb_w` (finite; the last one Z only).

Setup: 13 TeV, PDF 325100, fixed scale 91.1876 GeV, G_mu scheme, inputs of `base.input`.
Cuts (analysis `V_lljj_4q` here, `pwhg_analysis-NSCjj.f` in POWHEG): the two final partons with pT > 30 each and
Delta R > 0.4, no rapidity cut; lepton cuts as in `../llj`.

## Run

    # mixed_ll, built with analysis=V_lljj_4q, from f90/
    tests/powheg/lljj/lljjrun.sh wp out_wp.txt tests/powheg/lljj/sectors_w.txt 2000000 8000000      # also wm; z with sectors_z.txt
    # POWHEG (build and keys: README in NSC-lnujjMIX), about 10 minutes per seed
    ../pwgrun.sh p_1 <path>/NSC-lnujjMIX 24 chooseDY=1 ncall1=500000 itmx1=3 ncall2=1000000 itmx2=5
    ../pwgcomb.py p_1 rate_os; ../pwgcomb.py p_1 rate_ss

The counterterms of the singular sectors do not pass the cuts (they have one parton less).

## Results (2026-11-14)

| | initial state | POWHEG [fb] | mixed_ll [fb] | difference [fb] |
|---|---|---|---|---|
| W+ | opposite sign | -0.5683 +- 0.0008 | -0.5690 +- 0.0010 | -0.0008 +- 0.0013 |
| W+ | same sign | 0.4628 +- 0.0042 | 0.4598 +- 0.0013 | -0.0031 +- 0.0044 |
| W- | opposite sign | -0.3864 +- 0.0008 | -0.3862 +- 0.0009 | 0.0002 +- 0.0012 |
| W- | same sign | 0.2644 +- 0.0035 | 0.2614 +- 0.0006 | -0.0030 +- 0.0036 |
| Z | opposite sign | -0.4496 +- 0.0006 | -0.4488 +- 0.0007 | 0.0007 +- 0.0009 |
| Z | same sign | 0.3080 +- 0.0028 | 0.3086 +- 0.0008 | 0.0006 +- 0.0030 |

Opposite sign: `ns_qqb` + `ns_qqb_w` (+ `ns_qqpb_w` for Z), `rate_os` in POWHEG; same sign: `ns_qq` + `ns_qqp_w`, `rate_ss`.
The two have opposite signs and largely cancel in the sum, so they are compared separately.
POWHEG: plain mean over 44 (W+) and 24 (W-, Z) seeds. The same-sign rate has a large-mjj tail with a few large weights;
the error-weighted mean (`../pwgw.py`) is 2% lower and biased. Lowering the generation cut `genptj` from 10 to 5 changes nothing.
mixed_ll: 300k / 1M points (singular sectors), 2M / 8M (finite channels).

Found with this comparison: in the Z finite channels the OpenLoops processes were registered with u and d exchanged
(PDG 2 used for d), so each amplitude was folded with the luminosity of the other isospin. Fixed in `mod_ol_interface.f90`;
before the fix Z opposite sign was -0.4513 +- 0.0007 (`ns_qqb_w` -0.1592, `ns_qqpb_w` -0.1114), same sign 0.3113 (`ns_qqp_w` 0.2091).
The shifts of the three channels are 3 to 6 sigma; the comparison of the sum alone would not have been decisive (1.9 sigma).

mixed_ll sectors:

| ch sec | W+ [fb] | W- [fb] | Z [fb] |
|---|---|---|---|
| ns_qq rr_5161a | 0.0237 +- 0.0001 | 0.0243 +- 0.0001 | 0.0177 +- 0.0001 |
| ns_qq rr_5161c | 0.0580 +- 0.0004 | 0.0562 +- 0.0004 | 0.0331 +- 0.0003 |
| ns_qq rr_5262a | 0.0236 +- 0.0001 | 0.0242 +- 0.0002 | 0.0177 +- 0.0002 |
| ns_qq rr_5262c | 0.0568 +- 0.0003 | 0.0568 +- 0.0003 | 0.0337 +- 0.0004 |
| ns_qqb rr_5161a | -0.0665 +- 0.0002 | -0.0354 +- 0.0001 | -0.0308 +- 0.0002 |
| ns_qqb rr_5161c | -0.1389 +- 0.0005 | -0.0707 +- 0.0003 | -0.0596 +- 0.0004 |
| ns_qqb rr_5262a | -0.0669 +- 0.0002 | -0.0357 +- 0.0001 | -0.0310 +- 0.0001 |
| ns_qqb rr_5262c | -0.1403 +- 0.0005 | -0.0709 +- 0.0003 | -0.0594 +- 0.0003 |
| ns_qqb_w rr | -0.1565 +- 0.0007 | -0.1736 +- 0.0008 | -0.1604 +- 0.0001 |
| ns_qqp_w rr | 0.2977 +- 0.0012 | 0.0999 +- 0.0003 | 0.2063 +- 0.0006 |
| ns_qqpb_w rr | - | - | -0.1076 +- 0.0004 |
