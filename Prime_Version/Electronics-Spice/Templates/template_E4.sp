* =====================================================================
* Template for E4: 6T FinFET SRAM (22nm)
* =====================================================================
* FINAL VERSION:
*   1. .ic command placed in main netlist (outside .control) with {VDD}
*   2. Transient with 'uic' to force valid state (Q=VDD, Qb=0)
*   3. Leakage measured via 'let' vector + 'meas ... AVG' from 4n to 5n
*   4. Sentinel marker 'E4_SIMULATION_DONE' for completion verification
*   5. Comments use separate lines (no inline '*' comments)
* =====================================================================
.include "<PARAM_FILE>"
.include "<MODELS_DIR>/ptm22hp.lib"
.param pi=3.14159265359
.options reltol=1e-3 gmin=1e-12 itl1=100 itl2=50 itl4=50
.options method=gear cshunt=1e-15

* Physical fin width: 40nm per fin
.param W_single_fin = 40n

* Compute effective widths
.param W_pd = {Wfin_pd * W_single_fin}
.param W_pu = {Wfin_pu * W_single_fin}
.param W_ax = {Wfin_ax * W_single_fin}

* ---- Power supply and biases (Standby mode: WL=0) ----
VDD_SRC VDD 0 DC {VDD}

* Wordline disabled for standby leakage
VWL WL 0 DC 0

VBL BL 0 DC {VDD}
VBLB BLB 0 DC {VDD}

* ---- Transistors ----
Mpd1 Q Qb 0 0 nmos W={W_pd} L={Lfin_pd}
Mpu1 Q Qb VDD VDD pmos W={W_pu} L={Lfin_pu}
Mpd2 Qb Q 0 0 nmos W={W_pd} L={Lfin_pd}
Mpu2 Qb Q VDD VDD pmos W={W_pu} L={Lfin_pu}
Max1 BL WL Q 0 nmos W={W_ax} L={Lfin_ax}
Max2 BLB WL Qb 0 nmos W={W_ax} L={Lfin_ax}

* ---- Initial condition to force valid state (Q=VDD, Qb=0) ----
.ic V(Q)={VDD} V(Qb)=0

.control
  op
  echo ---- DC OPERATING POINT (STANDBY) ----
  print v(Q) v(Qb) v(BL) v(BLB) v(WL) v(VDD)
  echo -------------------------------------

  * ---- Transient with UIC to settle the state ----
  tran 0.05n 5n uic

  * ---- Leakage measurement: average supply current over last 1ns ----
  let i_leak_inst = abs(i(VDD_SRC))
  meas tran avg_leak AVG i_leak_inst FROM=4n TO=5n
  let I_leak = avg_leak
  echo I_leak_A = $&I_leak

  * ---- Sentinel marker ----
  echo E4_SIMULATION_DONE
.endc
.end