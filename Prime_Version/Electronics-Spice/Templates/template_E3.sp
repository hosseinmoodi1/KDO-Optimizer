* =====================================================================
* Template for E3: LC-VCO (65nm) - OPTIMIZED FOR SPEED
* =====================================================================
* OPTIMIZATIONS:
*   1. Reduced tstop from 500n to 300n (steady-state achieved earlier)
*   2. Increased maxstep from 1e-11 to 5e-11 (fewer simulation points)
*   3. Removed waveform saving (no .save or .write)
*   4. Reduced measurement window to 200n-300n
* =====================================================================
.include "<PARAM_FILE>"
.include "<MODELS_DIR>/ptm65.lib"
.param pi=3.14159265359

* ---- Speed-optimized options ----
.options reltol=1e-3 abstol=1e-12 vntol=1e-6
.options method=gear cshunt=2e-15 numdgt=5
.options maxstep=5e-11          * Increased from 1e-11
.options noacct                * Disable accounting (minor speedup)

* ---- Power Supply and Bias ----
VDD_SRC VDD 0 DC 1.2
I_tail tail 0 DC {Ibias}
Ctail tail 0 10p

* ---- Tank Inductor ----
Ltank     Voutp n_tank {L}
Rser_tank n_tank Voutn 1.5

* ---- Cross-Coupled Transistors ----
Mp1 Voutp Voutn VDD VDD pmos W={Wp} L=0.065u
Mp2 Voutn Voutp VDD VDD pmos W={Wp} L=0.065u
Mn1 Voutp Voutn tail 0 nmos W={Wn} L=0.065u
Mn2 Voutn Voutp tail 0 nmos W={Wn} L=0.065u

* ---- Symmetric Tank Capacitors ----
Ctank  Voutp Voutn {Cvar}
Cfixp  Voutp 0 {Cfix}
Cfixn  Voutn 0 {Cfix}
Cbufp  Voutp 0 {Cbuf}
Cbufn  Voutn 0 {Cbuf}

* ---- Startup Kick ----
Isync Voutp Voutn PULSE(0 1m 0 1n 1n 2n 10u)

* ---- Analytical Parameters (Not Connected) ----
* Rbias, Vctrl_min, Vctrl_max, K_vco

.control
  * ---- Disable output to reduce disk I/O ----
  set noaskquit
  set filetype=ascii
  * ---- Suppress echo of commands ----
  set quiet

  op
  echo ---- DC OPERATING POINT ----
  print v(Voutp) v(Voutn) v(tail)
  echo -----------------------------

  * ---- Reduced transient (300ns instead of 500ns) ----
  * Startup settles within ~100ns, measurement from 200ns to 300ns
  tran 0.05n 300n

  let p_inst = -v(VDD) * i(VDD_SRC)
  meas tran avg_power AVG p_inst FROM=200n TO=300n
  let Power_W = avg_power
  echo Power_W = $&Power_W

  let diff_out = v(Voutp) - v(Voutn)
  meas tran tperiod TRIG diff_out VAL=0 TD=200n RISE=1 TARG diff_out VAL=0 TD=200n RISE=2
  echo tperiod = $&tperiod

  meas tran diff_max MAX diff_out FROM=200n TO=300n
  meas tran diff_min MIN diff_out FROM=200n TO=300n
  let vpp_out = diff_max - diff_min
  echo vpp_out = $&vpp_out

  meas tran voutp_max MAX v(Voutp) FROM=200n TO=300n
  meas tran voutp_min MIN v(Voutp) FROM=200n TO=300n
  meas tran voutn_max MAX v(Voutn) FROM=200n TO=300n
  meas tran voutn_min MIN v(Voutn) FROM=200n TO=300n

  echo voutp_max = $&voutp_max
  echo voutp_min = $&voutp_min
  echo voutn_max = $&voutn_max
  echo voutn_min = $&voutn_min

  echo E3_SIMULATION_DONE
.endc
.end