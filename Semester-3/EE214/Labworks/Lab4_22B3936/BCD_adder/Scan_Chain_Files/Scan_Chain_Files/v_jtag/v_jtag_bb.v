
module v_jtag (
	virtual_jtag_tdi,
	virtual_jtag_tdo,
	virtual_jtag_ir_in,
	virtual_jtag_ir_out,
	virtual_jtag_virtual_state_cdr,
	virtual_jtag_virtual_state_sdr,
	virtual_jtag_virtual_state_e1dr,
	virtual_jtag_virtual_state_pdr,
	virtual_jtag_virtual_state_e2dr,
	virtual_jtag_virtual_state_udr,
	virtual_jtag_virtual_state_cir,
	virtual_jtag_virtual_state_uir,
	virtual_jtag_tms,
	virtual_jtag_jtag_state_tlr,
	virtual_jtag_jtag_state_rti,
	virtual_jtag_jtag_state_sdrs,
	virtual_jtag_jtag_state_cdr,
	virtual_jtag_jtag_state_sdr,
	virtual_jtag_jtag_state_e1dr,
	virtual_jtag_jtag_state_pdr,
	virtual_jtag_jtag_state_e2dr,
	virtual_jtag_jtag_state_udr,
	virtual_jtag_jtag_state_sirs,
	virtual_jtag_jtag_state_cir,
	virtual_jtag_jtag_state_sir,
	virtual_jtag_jtag_state_e1ir,
	virtual_jtag_jtag_state_pir,
	virtual_jtag_jtag_state_e2ir,
	virtual_jtag_jtag_state_uir,
	tck_clk);	

	output		virtual_jtag_tdi;
	input		virtual_jtag_tdo;
	output	[0:0]	virtual_jtag_ir_in;
	input	[0:0]	virtual_jtag_ir_out;
	output		virtual_jtag_virtual_state_cdr;
	output		virtual_jtag_virtual_state_sdr;
	output		virtual_jtag_virtual_state_e1dr;
	output		virtual_jtag_virtual_state_pdr;
	output		virtual_jtag_virtual_state_e2dr;
	output		virtual_jtag_virtual_state_udr;
	output		virtual_jtag_virtual_state_cir;
	output		virtual_jtag_virtual_state_uir;
	output		virtual_jtag_tms;
	output		virtual_jtag_jtag_state_tlr;
	output		virtual_jtag_jtag_state_rti;
	output		virtual_jtag_jtag_state_sdrs;
	output		virtual_jtag_jtag_state_cdr;
	output		virtual_jtag_jtag_state_sdr;
	output		virtual_jtag_jtag_state_e1dr;
	output		virtual_jtag_jtag_state_pdr;
	output		virtual_jtag_jtag_state_e2dr;
	output		virtual_jtag_jtag_state_udr;
	output		virtual_jtag_jtag_state_sirs;
	output		virtual_jtag_jtag_state_cir;
	output		virtual_jtag_jtag_state_sir;
	output		virtual_jtag_jtag_state_e1ir;
	output		virtual_jtag_jtag_state_pir;
	output		virtual_jtag_jtag_state_e2ir;
	output		virtual_jtag_jtag_state_uir;
	output		tck_clk;
endmodule
