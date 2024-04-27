transcript on
if {[file exists rtl_work]} {
	vdel -lib rtl_work -all
}
vlib rtl_work
vmap work rtl_work

vcom -93 -work work {C:/Users/ASUS/OneDrive - Indian Institute of Technology Bombay/ACADS/SEMESTER-3/EE-214/Practice_Problems/Level-III/STUDENTS_FSM_AmanR/STUDENTS_FSM/Gates.vhdl}
vcom -93 -work work {C:/Users/ASUS/OneDrive - Indian Institute of Technology Bombay/ACADS/SEMESTER-3/EE-214/Practice_Problems/Level-III/STUDENTS_FSM_AmanR/STUDENTS_FSM/AND_9.vhd}
vcom -93 -work work {C:/Users/ASUS/OneDrive - Indian Institute of Technology Bombay/ACADS/SEMESTER-3/EE-214/Practice_Problems/Level-III/STUDENTS_FSM_AmanR/STUDENTS_FSM/J_KFF.vhd}
vcom -93 -work work {C:/Users/ASUS/OneDrive - Indian Institute of Technology Bombay/ACADS/SEMESTER-3/EE-214/Practice_Problems/Level-III/STUDENTS_FSM_AmanR/STUDENTS_FSM/DFF_JK.vhd}
vcom -93 -work work {C:/Users/ASUS/OneDrive - Indian Institute of Technology Bombay/ACADS/SEMESTER-3/EE-214/Practice_Problems/Level-III/STUDENTS_FSM_AmanR/STUDENTS_FSM/STUDENTS_FSM.vhd}

vcom -93 -work work {C:/Users/ASUS/OneDrive - Indian Institute of Technology Bombay/ACADS/SEMESTER-3/EE-214/Practice_Problems/Level-III/STUDENTS_FSM_AmanR/STUDENTS_FSM/Testbench.vhdl}

vsim -t 1ps -L altera -L lpm -L sgate -L altera_mf -L altera_lnsim -L fiftyfivenm -L rtl_work -L work -voptargs="+acc"  Testbench

add wave *
view structure
view signals
run -all
