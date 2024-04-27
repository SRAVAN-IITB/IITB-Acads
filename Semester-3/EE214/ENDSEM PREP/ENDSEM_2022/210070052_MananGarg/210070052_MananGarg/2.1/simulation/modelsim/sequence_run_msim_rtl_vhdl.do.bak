transcript on
if {[file exists rtl_work]} {
	vdel -lib rtl_work -all
}
vlib rtl_work
vmap work rtl_work

vcom -93 -work work {E:/Sem-3/EE214/Finals/210070052_MananGarg/2.1/Sequence.vhdl}
vcom -93 -work work {E:/Sem-3/EE214/Finals/210070052_MananGarg/2.1/DUT.vhdl}

vcom -93 -work work {E:/Sem-3/EE214/Finals/210070052_MananGarg/2.1/Testbench.vhdl}

vsim -t 1ps -L altera -L lpm -L sgate -L altera_mf -L altera_lnsim -L fiftyfivenm -L rtl_work -L work -voptargs="+acc"  testbench

add wave *
view structure
view signals
run -all
