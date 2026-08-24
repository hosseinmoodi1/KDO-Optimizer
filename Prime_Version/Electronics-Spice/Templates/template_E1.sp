* =====================================================================
* Template for E1: Two-Stage Miller-Compensated Op-Amp (180nm)
* =====================================================================
* FINAL VERSION: Robust PM measurement using GBW scalar with at=$&GBW.
* =====================================================================
.include "<PARAM_FILE>"
.include "<MODELS_DIR>/ptm180.lib"

.options reltol=1e-2 gmin=1e-12 itl1=200 itl2=100 itl4=100
.options method=trap cshunt=1e-15 numdgt=5

VDD_SRC VDD 0 DC 1.8
VINP VINP 0 DC 0.9 AC 1

R_fb Vout VINN 1G
C_ac VINN 0 1G

Itail tail 0 DC {Ibias}
I_bias_ref Vb_ref 0 DC {Ibias/2}
M8 Vb_ref Vb_ref VDD VDD pmos W={W8} L={L8}

M1 n1 VINP tail 0 nmos W={W1} L={L1}
M2 n2 VINN tail 0 nmos W={W2} L={L2}

M3 n1 n1 VDD VDD pmos W={W3} L={L3}
M4 n2 n1 VDD VDD pmos W={W4} L={L4}

M6 Vout n2 0 0 nmos W={W6} L={L6}
M7 Vout Vb_ref VDD VDD pmos W={W7} L={L7}

CCOMP n2 Vout {Cc}
CL Vout 0 2p

.control
  set noaskquit
  op

  echo ---- DC NODE VOLTAGES ----
  print v(Vout) v(n1) v(n2) v(tail) v(Vb_ref)
  echo --------------------------

  let Vout_dc = v(Vout)
  let Vn1_dc = v(n1)
  let Vn2_dc = v(n2)
  let Vtail_dc = v(tail)
  let Vb_ref_dc = v(Vb_ref)
  echo Vout_dc = $&Vout_dc
  echo Vn1_dc = $&Vn1_dc
  echo Vn2_dc = $&Vn2_dc
  echo Vtail_dc = $&Vtail_dc
  echo Vb_ref_dc = $&Vb_ref_dc

  let Power_W = -v(VDD)*i(VDD_SRC)
  echo Power_W = $&Power_W

  ac dec 200 1 100000meg

  let gain = v(Vout) / v(VINP)
  let gain_db = vdb(gain)

  meas ac Gain_low_db find gain_db at=1
  echo Gain_low_db = $&Gain_low_db

  meas ac Gain_max_db MAX gain_db
  echo Gain_max_db = $&Gain_max_db

  * 1. Measure GBW when gain crosses 0 dB (falling edge)
  meas ac gbw when gain_db=0 fall=1

  if $?gbw
    let GBW = gbw
    echo GBW = $&GBW
  else
    if Gain_low_db > 0
      let target_gain = Gain_low_db - 3
      meas ac gbw_3db when gain_db=target_gain fall=1
      if $?gbw_3db
        let GBW = gbw_3db
        echo GBW = $&GBW (approx)
      else
        let GBW = 0
        echo GBW = 0
      end
    else
      let GBW = 0
      echo GBW = 0
    end
  end

  * 2. Measure PM robustly using the found GBW scalar
  * This bypasses the WHEN bug in newer ngspice versions
  let PM = -180
  if Gain_low_db > 0
    if $?gbw
      * $&GBW passes the scalar frequency directly to the at= parameter
      meas ac phase_at_gbw FIND ph(gain) at=$&GBW
      if $?phase_at_gbw
        let PM = 180 - abs(phase_at_gbw)
      end
    end
  end
  echo PM = $&PM

  let phase_1hz = ph(gain)[1]
  let phase_1mhz = ph(gain)[1e6]
  echo Phase_1Hz = $&phase_1hz
  echo Phase_1MHz = $&phase_1mhz
.endc
.end