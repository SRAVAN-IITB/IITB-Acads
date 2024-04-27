transcript on
if {[file exists rtl_work]} {
	vdel -lib rtl_work -all
}
vlib rtl_work
vmap work rtl_work

vcom -93 -work work {C:/Users/ASUS/OneDrive - Indian Institute of Technology Bombay/ACADS/SEMESTER-3/EE-214/ENDSEM PRACTICE/BCDtoSeven-Segment_Display/BCDtoSevenSegment_Display_bhv_usingCASE.vhdl}

vcom -93 -work work {C:/Users/ASUS/OneDrive - Indian Institute of Technology Bombay/ACADS/SEMESTER-3/EE-214/ENDSEM PRACTICE/BCDtoSeven-Segment_Display/BCDtoSevenSegment_Display_bhv_usingCASE.vhdl}
vcom -93 -work work {C:/Users/ASUS/OneDrive - Indian Institute of Technology Bombay/ACADS/SEMESTER-3/EE-214/ENDSEM PRACTICE/BCDtoSeven-Segment_Display/tb_bcd_7seg.vhdl}

vsim -t 1ps -L altera -L lpm -L sgate -L altera_mf -L altera_lnsim -L fiftyfivenm -L rtl_work -L work -voptargs="+acc"  tb_bcd_7seg

add wave *
view structure
view signals
run -all
