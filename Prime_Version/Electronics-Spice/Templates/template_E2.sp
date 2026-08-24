* =====================================================================
* Template for E2: RF Cascode LNA (65nm)
* =====================================================================
* FIXED:
*   1. Lg2 connected to ground (0)
*   2. Cascode gate biased via resistor from diode network
*   3. Decoupling capacitor for cascode gate
*   4. DC operating point printed for diagnostics
* =====================================================================
.include "<PARAM_FILE>"
.include "<MODEL_FILE>"
.param pi=3.14159265359
.options reltol=1e-2 gmin=1e-12 itl1=100 itl2=50 itl4=50
.options method=trap cshunt=1e-15 numdgt=5

VDD_SRC Vdd 0 DC {Vdd}
Vin node_rf_in 0 DC 0 AC 1
Rsin node_rf_in node_gate_in 50
Lg1 node_gate_in node_gate {Lg1}
Lg2 node_source 0 {Lg2}
Rg node_gate node_bias_d 10k
Ld Vdd node_drain {Ld}

M1 node_cascode_src node_gate node_source 0 nmos W={Wm1} L=0.065u
M2 node_drain node_casc node_cascode_src 0 nmos W={Wm2} L=0.065u

Cgd1 node_gate node_cascode_src {Cgd1}
Cgd2 node_gate node_drain {Cgd2}
Cgd3 node_casc node_drain {Cgd3}

M3 node_bias_d node_bias_d 0 0 nmos W={Wm3} L=0.065u
M4 node_bias_d node_bias_d 0 0 nmos W={Wm4} L=0.065u

Ibias Vdd node_bias_d DC {Id}

Rbias node_bias_d node_casc 10k
Ccasc node_casc 0 10p

Cout node_drain node_out 1p
Rload node_out 0 50

.control
  op
  echo ---- DC OPERATING POINT ----
  print v(node_bias_d) v(node_casc) v(node_gate) v(node_source) v(node_cascode_src) v(node_drain)
  echo -----------------------------

  let Power_W = -v(Vdd)*i(VDD_SRC)
  echo Power_W = $&Power_W

  ac dec 200 1e9 10e9
  meas ac Gain_db max vdb(node_out) from=2.3e9 to=2.5e9
  echo Gain_db = $&Gain_db
.endc
.end