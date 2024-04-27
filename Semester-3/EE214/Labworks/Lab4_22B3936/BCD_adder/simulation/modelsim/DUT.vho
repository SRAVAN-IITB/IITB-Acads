-- Copyright (C) 2018  Intel Corporation. All rights reserved.
-- Your use of Intel Corporation's design tools, logic functions 
-- and other software and tools, and its AMPP partner logic 
-- functions, and any output files from any of the foregoing 
-- (including device programming or simulation files), and any 
-- associated documentation or information are expressly subject 
-- to the terms and conditions of the Intel Program License 
-- Subscription Agreement, the Intel Quartus Prime License Agreement,
-- the Intel FPGA IP License Agreement, or other applicable license
-- agreement, including, without limitation, that your use is for
-- the sole purpose of programming logic devices manufactured by
-- Intel and sold by Intel or its authorized distributors.  Please
-- refer to the applicable agreement for further details.

-- VENDOR "Altera"
-- PROGRAM "Quartus Prime"
-- VERSION "Version 18.1.0 Build 625 09/12/2018 SJ Lite Edition"

-- DATE "09/01/2023 17:16:32"

-- 
-- Device: Altera 10M25SAE144C8G Package EQFP144
-- 

-- 
-- This VHDL file should be used for ModelSim-Altera (VHDL) only
-- 

LIBRARY FIFTYFIVENM;
LIBRARY IEEE;
USE FIFTYFIVENM.FIFTYFIVENM_COMPONENTS.ALL;
USE IEEE.STD_LOGIC_1164.ALL;

ENTITY 	hard_block IS
    PORT (
	devoe : IN std_logic;
	devclrn : IN std_logic;
	devpor : IN std_logic
	);
END hard_block;

-- Design Ports Information
-- ~ALTERA_TMS~	=>  Location: PIN_16,	 I/O Standard: 2.5 V Schmitt Trigger,	 Current Strength: Default
-- ~ALTERA_TCK~	=>  Location: PIN_18,	 I/O Standard: 2.5 V Schmitt Trigger,	 Current Strength: Default
-- ~ALTERA_TDI~	=>  Location: PIN_19,	 I/O Standard: 2.5 V Schmitt Trigger,	 Current Strength: Default
-- ~ALTERA_TDO~	=>  Location: PIN_20,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- ~ALTERA_CONFIG_SEL~	=>  Location: PIN_126,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- ~ALTERA_nCONFIG~	=>  Location: PIN_129,	 I/O Standard: 2.5 V Schmitt Trigger,	 Current Strength: Default
-- ~ALTERA_nSTATUS~	=>  Location: PIN_136,	 I/O Standard: 2.5 V Schmitt Trigger,	 Current Strength: Default
-- ~ALTERA_CONF_DONE~	=>  Location: PIN_138,	 I/O Standard: 2.5 V Schmitt Trigger,	 Current Strength: Default


ARCHITECTURE structure OF hard_block IS
SIGNAL gnd : std_logic := '0';
SIGNAL vcc : std_logic := '1';
SIGNAL unknown : std_logic := 'X';
SIGNAL ww_devoe : std_logic;
SIGNAL ww_devclrn : std_logic;
SIGNAL ww_devpor : std_logic;
SIGNAL \~ALTERA_TMS~~padout\ : std_logic;
SIGNAL \~ALTERA_TCK~~padout\ : std_logic;
SIGNAL \~ALTERA_TDI~~padout\ : std_logic;
SIGNAL \~ALTERA_CONFIG_SEL~~padout\ : std_logic;
SIGNAL \~ALTERA_nCONFIG~~padout\ : std_logic;
SIGNAL \~ALTERA_nSTATUS~~padout\ : std_logic;
SIGNAL \~ALTERA_CONF_DONE~~padout\ : std_logic;
SIGNAL \~ALTERA_TMS~~ibuf_o\ : std_logic;
SIGNAL \~ALTERA_TCK~~ibuf_o\ : std_logic;
SIGNAL \~ALTERA_TDI~~ibuf_o\ : std_logic;
SIGNAL \~ALTERA_CONFIG_SEL~~ibuf_o\ : std_logic;
SIGNAL \~ALTERA_nCONFIG~~ibuf_o\ : std_logic;
SIGNAL \~ALTERA_nSTATUS~~ibuf_o\ : std_logic;
SIGNAL \~ALTERA_CONF_DONE~~ibuf_o\ : std_logic;

BEGIN

ww_devoe <= devoe;
ww_devclrn <= devclrn;
ww_devpor <= devpor;
END structure;


LIBRARY FIFTYFIVENM;
LIBRARY IEEE;
USE FIFTYFIVENM.FIFTYFIVENM_COMPONENTS.ALL;
USE IEEE.STD_LOGIC_1164.ALL;

ENTITY 	DUT IS
    PORT (
	input_vector : IN std_logic_vector(7 DOWNTO 0);
	output_vector : OUT std_logic_vector(4 DOWNTO 0)
	);
END DUT;

-- Design Ports Information
-- output_vector[0]	=>  Location: PIN_124,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- output_vector[1]	=>  Location: PIN_135,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- output_vector[2]	=>  Location: PIN_22,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- output_vector[3]	=>  Location: PIN_127,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- output_vector[4]	=>  Location: PIN_25,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- input_vector[0]	=>  Location: PIN_131,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- input_vector[4]	=>  Location: PIN_132,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- input_vector[1]	=>  Location: PIN_123,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- input_vector[5]	=>  Location: PIN_140,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- input_vector[2]	=>  Location: PIN_134,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- input_vector[6]	=>  Location: PIN_21,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- input_vector[7]	=>  Location: PIN_141,	 I/O Standard: 2.5 V,	 Current Strength: Default
-- input_vector[3]	=>  Location: PIN_130,	 I/O Standard: 2.5 V,	 Current Strength: Default


ARCHITECTURE structure OF DUT IS
SIGNAL gnd : std_logic := '0';
SIGNAL vcc : std_logic := '1';
SIGNAL unknown : std_logic := 'X';
SIGNAL devoe : std_logic := '1';
SIGNAL devclrn : std_logic := '1';
SIGNAL devpor : std_logic := '1';
SIGNAL ww_devoe : std_logic;
SIGNAL ww_devclrn : std_logic;
SIGNAL ww_devpor : std_logic;
SIGNAL ww_input_vector : std_logic_vector(7 DOWNTO 0);
SIGNAL ww_output_vector : std_logic_vector(4 DOWNTO 0);
SIGNAL \~QUARTUS_CREATED_ADC1~_CHSEL_bus\ : std_logic_vector(4 DOWNTO 0);
SIGNAL \~QUARTUS_CREATED_ADC2~_CHSEL_bus\ : std_logic_vector(4 DOWNTO 0);
SIGNAL \~QUARTUS_CREATED_GND~I_combout\ : std_logic;
SIGNAL \~QUARTUS_CREATED_UNVM~~busy\ : std_logic;
SIGNAL \~QUARTUS_CREATED_ADC1~~eoc\ : std_logic;
SIGNAL \~QUARTUS_CREATED_ADC2~~eoc\ : std_logic;
SIGNAL \output_vector[0]~output_o\ : std_logic;
SIGNAL \output_vector[1]~output_o\ : std_logic;
SIGNAL \output_vector[2]~output_o\ : std_logic;
SIGNAL \output_vector[3]~output_o\ : std_logic;
SIGNAL \output_vector[4]~output_o\ : std_logic;
SIGNAL \input_vector[0]~input_o\ : std_logic;
SIGNAL \input_vector[4]~input_o\ : std_logic;
SIGNAL \add_instance|fourbit_adder_0|Half_adder_0|S~combout\ : std_logic;
SIGNAL \input_vector[1]~input_o\ : std_logic;
SIGNAL \input_vector[5]~input_o\ : std_logic;
SIGNAL \add_instance|fourbit_adder_0|Full_adder_1|XOR2|Y~0_combout\ : std_logic;
SIGNAL \input_vector[7]~input_o\ : std_logic;
SIGNAL \add_instance|fourbit_adder_0|Full_adder_1|OR2|Y~0_combout\ : std_logic;
SIGNAL \input_vector[6]~input_o\ : std_logic;
SIGNAL \input_vector[2]~input_o\ : std_logic;
SIGNAL \add_instance|fourbit_adder_0|Full_adder_2|OR2|Y~0_combout\ : std_logic;
SIGNAL \add_instance|fourbit_adder_1|Full_adder_3|XOR2|Y~0_combout\ : std_logic;
SIGNAL \input_vector[3]~input_o\ : std_logic;
SIGNAL \add_instance|OR2|Y~0_combout\ : std_logic;
SIGNAL \add_instance|fourbit_adder_1|Full_adder_1|XOR1|Y~combout\ : std_logic;
SIGNAL \add_instance|fourbit_adder_0|Full_adder_2|XOR2|Y~0_combout\ : std_logic;
SIGNAL \add_instance|fourbit_adder_1|Full_adder_2|XOR2|Y~0_combout\ : std_logic;
SIGNAL \add_instance|fourbit_adder_1|Full_adder_3|XOR2|Y~1_combout\ : std_logic;

COMPONENT hard_block
    PORT (
	devoe : IN std_logic;
	devclrn : IN std_logic;
	devpor : IN std_logic);
END COMPONENT;

BEGIN

ww_input_vector <= input_vector;
output_vector <= ww_output_vector;
ww_devoe <= devoe;
ww_devclrn <= devclrn;
ww_devpor <= devpor;

\~QUARTUS_CREATED_ADC1~_CHSEL_bus\ <= (\~QUARTUS_CREATED_GND~I_combout\ & \~QUARTUS_CREATED_GND~I_combout\ & \~QUARTUS_CREATED_GND~I_combout\ & \~QUARTUS_CREATED_GND~I_combout\ & \~QUARTUS_CREATED_GND~I_combout\);

\~QUARTUS_CREATED_ADC2~_CHSEL_bus\ <= (\~QUARTUS_CREATED_GND~I_combout\ & \~QUARTUS_CREATED_GND~I_combout\ & \~QUARTUS_CREATED_GND~I_combout\ & \~QUARTUS_CREATED_GND~I_combout\ & \~QUARTUS_CREATED_GND~I_combout\);
auto_generated_inst : hard_block
PORT MAP (
	devoe => ww_devoe,
	devclrn => ww_devclrn,
	devpor => ww_devpor);

-- Location: LCCOMB_X26_Y23_N24
\~QUARTUS_CREATED_GND~I\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \~QUARTUS_CREATED_GND~I_combout\ = GND

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0000000000000000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	combout => \~QUARTUS_CREATED_GND~I_combout\);

-- Location: IOOBUF_X14_Y21_N23
\output_vector[0]~output\ : fiftyfivenm_io_obuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	open_drain_output => "false")
-- pragma translate_on
PORT MAP (
	i => \add_instance|fourbit_adder_0|Half_adder_0|S~combout\,
	devoe => ww_devoe,
	o => \output_vector[0]~output_o\);

-- Location: IOOBUF_X10_Y21_N16
\output_vector[1]~output\ : fiftyfivenm_io_obuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	open_drain_output => "false")
-- pragma translate_on
PORT MAP (
	i => \add_instance|fourbit_adder_1|Full_adder_1|XOR1|Y~combout\,
	devoe => ww_devoe,
	o => \output_vector[1]~output_o\);

-- Location: IOOBUF_X0_Y16_N23
\output_vector[2]~output\ : fiftyfivenm_io_obuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	open_drain_output => "false")
-- pragma translate_on
PORT MAP (
	i => \add_instance|fourbit_adder_1|Full_adder_2|XOR2|Y~0_combout\,
	devoe => ww_devoe,
	o => \output_vector[2]~output_o\);

-- Location: IOOBUF_X14_Y21_N30
\output_vector[3]~output\ : fiftyfivenm_io_obuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	open_drain_output => "false")
-- pragma translate_on
PORT MAP (
	i => \add_instance|fourbit_adder_1|Full_adder_3|XOR2|Y~1_combout\,
	devoe => ww_devoe,
	o => \output_vector[3]~output_o\);

-- Location: IOOBUF_X0_Y15_N23
\output_vector[4]~output\ : fiftyfivenm_io_obuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	open_drain_output => "false")
-- pragma translate_on
PORT MAP (
	i => \add_instance|OR2|Y~0_combout\,
	devoe => ww_devoe,
	o => \output_vector[4]~output_o\);

-- Location: IOIBUF_X12_Y21_N8
\input_vector[0]~input\ : fiftyfivenm_io_ibuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	listen_to_nsleep_signal => "false",
	simulate_z_as => "z")
-- pragma translate_on
PORT MAP (
	i => ww_input_vector(0),
	o => \input_vector[0]~input_o\);

-- Location: IOIBUF_X12_Y21_N29
\input_vector[4]~input\ : fiftyfivenm_io_ibuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	listen_to_nsleep_signal => "false",
	simulate_z_as => "z")
-- pragma translate_on
PORT MAP (
	i => ww_input_vector(4),
	o => \input_vector[4]~input_o\);

-- Location: LCCOMB_X8_Y17_N0
\add_instance|fourbit_adder_0|Half_adder_0|S\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \add_instance|fourbit_adder_0|Half_adder_0|S~combout\ = \input_vector[0]~input_o\ $ (\input_vector[4]~input_o\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0011001111001100",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \input_vector[0]~input_o\,
	datad => \input_vector[4]~input_o\,
	combout => \add_instance|fourbit_adder_0|Half_adder_0|S~combout\);

-- Location: IOIBUF_X14_Y21_N1
\input_vector[1]~input\ : fiftyfivenm_io_ibuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	listen_to_nsleep_signal => "false",
	simulate_z_as => "z")
-- pragma translate_on
PORT MAP (
	i => ww_input_vector(1),
	o => \input_vector[1]~input_o\);

-- Location: IOIBUF_X8_Y21_N15
\input_vector[5]~input\ : fiftyfivenm_io_ibuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	listen_to_nsleep_signal => "false",
	simulate_z_as => "z")
-- pragma translate_on
PORT MAP (
	i => ww_input_vector(5),
	o => \input_vector[5]~input_o\);

-- Location: LCCOMB_X8_Y17_N10
\add_instance|fourbit_adder_0|Full_adder_1|XOR2|Y~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \add_instance|fourbit_adder_0|Full_adder_1|XOR2|Y~0_combout\ = \input_vector[1]~input_o\ $ (\input_vector[5]~input_o\ $ (((\input_vector[0]~input_o\ & \input_vector[4]~input_o\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1001011001011010",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \input_vector[1]~input_o\,
	datab => \input_vector[0]~input_o\,
	datac => \input_vector[5]~input_o\,
	datad => \input_vector[4]~input_o\,
	combout => \add_instance|fourbit_adder_0|Full_adder_1|XOR2|Y~0_combout\);

-- Location: IOIBUF_X8_Y21_N22
\input_vector[7]~input\ : fiftyfivenm_io_ibuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	listen_to_nsleep_signal => "false",
	simulate_z_as => "z")
-- pragma translate_on
PORT MAP (
	i => ww_input_vector(7),
	o => \input_vector[7]~input_o\);

-- Location: LCCOMB_X8_Y17_N28
\add_instance|fourbit_adder_0|Full_adder_1|OR2|Y~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \add_instance|fourbit_adder_0|Full_adder_1|OR2|Y~0_combout\ = (\input_vector[1]~input_o\ & ((\input_vector[5]~input_o\) # ((\input_vector[0]~input_o\ & \input_vector[4]~input_o\)))) # (!\input_vector[1]~input_o\ & (\input_vector[0]~input_o\ & 
-- (\input_vector[5]~input_o\ & \input_vector[4]~input_o\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1110100010100000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \input_vector[1]~input_o\,
	datab => \input_vector[0]~input_o\,
	datac => \input_vector[5]~input_o\,
	datad => \input_vector[4]~input_o\,
	combout => \add_instance|fourbit_adder_0|Full_adder_1|OR2|Y~0_combout\);

-- Location: IOIBUF_X0_Y16_N15
\input_vector[6]~input\ : fiftyfivenm_io_ibuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	listen_to_nsleep_signal => "false",
	simulate_z_as => "z")
-- pragma translate_on
PORT MAP (
	i => ww_input_vector(6),
	o => \input_vector[6]~input_o\);

-- Location: IOIBUF_X10_Y21_N1
\input_vector[2]~input\ : fiftyfivenm_io_ibuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	listen_to_nsleep_signal => "false",
	simulate_z_as => "z")
-- pragma translate_on
PORT MAP (
	i => ww_input_vector(2),
	o => \input_vector[2]~input_o\);

-- Location: LCCOMB_X8_Y17_N24
\add_instance|fourbit_adder_0|Full_adder_2|OR2|Y~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \add_instance|fourbit_adder_0|Full_adder_2|OR2|Y~0_combout\ = (\add_instance|fourbit_adder_0|Full_adder_1|OR2|Y~0_combout\ & ((\input_vector[6]~input_o\) # (\input_vector[2]~input_o\))) # (!\add_instance|fourbit_adder_0|Full_adder_1|OR2|Y~0_combout\ & 
-- (\input_vector[6]~input_o\ & \input_vector[2]~input_o\))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111110011000000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \add_instance|fourbit_adder_0|Full_adder_1|OR2|Y~0_combout\,
	datac => \input_vector[6]~input_o\,
	datad => \input_vector[2]~input_o\,
	combout => \add_instance|fourbit_adder_0|Full_adder_2|OR2|Y~0_combout\);

-- Location: LCCOMB_X8_Y17_N22
\add_instance|fourbit_adder_1|Full_adder_3|XOR2|Y~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \add_instance|fourbit_adder_1|Full_adder_3|XOR2|Y~0_combout\ = (\add_instance|fourbit_adder_0|Full_adder_1|XOR2|Y~0_combout\) # (\add_instance|fourbit_adder_0|Full_adder_1|OR2|Y~0_combout\ $ (\input_vector[6]~input_o\ $ (\input_vector[2]~input_o\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1110101110111110",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \add_instance|fourbit_adder_0|Full_adder_1|XOR2|Y~0_combout\,
	datab => \add_instance|fourbit_adder_0|Full_adder_1|OR2|Y~0_combout\,
	datac => \input_vector[6]~input_o\,
	datad => \input_vector[2]~input_o\,
	combout => \add_instance|fourbit_adder_1|Full_adder_3|XOR2|Y~0_combout\);

-- Location: IOIBUF_X12_Y21_N1
\input_vector[3]~input\ : fiftyfivenm_io_ibuf
-- pragma translate_off
GENERIC MAP (
	bus_hold => "false",
	listen_to_nsleep_signal => "false",
	simulate_z_as => "z")
-- pragma translate_on
PORT MAP (
	i => ww_input_vector(3),
	o => \input_vector[3]~input_o\);

-- Location: LCCOMB_X8_Y17_N26
\add_instance|OR2|Y~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \add_instance|OR2|Y~0_combout\ = (\input_vector[7]~input_o\ & ((\add_instance|fourbit_adder_0|Full_adder_2|OR2|Y~0_combout\) # ((\add_instance|fourbit_adder_1|Full_adder_3|XOR2|Y~0_combout\) # (\input_vector[3]~input_o\)))) # (!\input_vector[7]~input_o\ & 
-- ((\add_instance|fourbit_adder_0|Full_adder_2|OR2|Y~0_combout\ & ((\add_instance|fourbit_adder_1|Full_adder_3|XOR2|Y~0_combout\) # (\input_vector[3]~input_o\))) # (!\add_instance|fourbit_adder_0|Full_adder_2|OR2|Y~0_combout\ & 
-- (\add_instance|fourbit_adder_1|Full_adder_3|XOR2|Y~0_combout\ & \input_vector[3]~input_o\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1111111011101000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \input_vector[7]~input_o\,
	datab => \add_instance|fourbit_adder_0|Full_adder_2|OR2|Y~0_combout\,
	datac => \add_instance|fourbit_adder_1|Full_adder_3|XOR2|Y~0_combout\,
	datad => \input_vector[3]~input_o\,
	combout => \add_instance|OR2|Y~0_combout\);

-- Location: LCCOMB_X8_Y17_N12
\add_instance|fourbit_adder_1|Full_adder_1|XOR1|Y\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \add_instance|fourbit_adder_1|Full_adder_1|XOR1|Y~combout\ = \add_instance|fourbit_adder_0|Full_adder_1|XOR2|Y~0_combout\ $ (\add_instance|OR2|Y~0_combout\)

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0101101001011010",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \add_instance|fourbit_adder_0|Full_adder_1|XOR2|Y~0_combout\,
	datac => \add_instance|OR2|Y~0_combout\,
	combout => \add_instance|fourbit_adder_1|Full_adder_1|XOR1|Y~combout\);

-- Location: LCCOMB_X8_Y17_N6
\add_instance|fourbit_adder_0|Full_adder_2|XOR2|Y~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \add_instance|fourbit_adder_0|Full_adder_2|XOR2|Y~0_combout\ = \add_instance|fourbit_adder_0|Full_adder_1|OR2|Y~0_combout\ $ (\input_vector[6]~input_o\ $ (\input_vector[2]~input_o\))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1100001100111100",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	datab => \add_instance|fourbit_adder_0|Full_adder_1|OR2|Y~0_combout\,
	datac => \input_vector[6]~input_o\,
	datad => \input_vector[2]~input_o\,
	combout => \add_instance|fourbit_adder_0|Full_adder_2|XOR2|Y~0_combout\);

-- Location: LCCOMB_X8_Y17_N16
\add_instance|fourbit_adder_1|Full_adder_2|XOR2|Y~0\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \add_instance|fourbit_adder_1|Full_adder_2|XOR2|Y~0_combout\ = \add_instance|fourbit_adder_0|Full_adder_2|XOR2|Y~0_combout\ $ (((!\add_instance|fourbit_adder_0|Full_adder_1|XOR2|Y~0_combout\ & \add_instance|OR2|Y~0_combout\)))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "1010111101010000",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \add_instance|fourbit_adder_0|Full_adder_1|XOR2|Y~0_combout\,
	datac => \add_instance|OR2|Y~0_combout\,
	datad => \add_instance|fourbit_adder_0|Full_adder_2|XOR2|Y~0_combout\,
	combout => \add_instance|fourbit_adder_1|Full_adder_2|XOR2|Y~0_combout\);

-- Location: LCCOMB_X8_Y17_N2
\add_instance|fourbit_adder_1|Full_adder_3|XOR2|Y~1\ : fiftyfivenm_lcell_comb
-- Equation(s):
-- \add_instance|fourbit_adder_1|Full_adder_3|XOR2|Y~1_combout\ = (\input_vector[7]~input_o\ & (\add_instance|fourbit_adder_0|Full_adder_2|OR2|Y~0_combout\ $ (\add_instance|fourbit_adder_1|Full_adder_3|XOR2|Y~0_combout\ $ (!\input_vector[3]~input_o\)))) # 
-- (!\input_vector[7]~input_o\ & ((\add_instance|fourbit_adder_0|Full_adder_2|OR2|Y~0_combout\ & (\add_instance|fourbit_adder_1|Full_adder_3|XOR2|Y~0_combout\ $ (!\input_vector[3]~input_o\))) # (!\add_instance|fourbit_adder_0|Full_adder_2|OR2|Y~0_combout\ & 
-- (!\add_instance|fourbit_adder_1|Full_adder_3|XOR2|Y~0_combout\ & \input_vector[3]~input_o\))))

-- pragma translate_off
GENERIC MAP (
	lut_mask => "0110100110000110",
	sum_lutc_input => "datac")
-- pragma translate_on
PORT MAP (
	dataa => \input_vector[7]~input_o\,
	datab => \add_instance|fourbit_adder_0|Full_adder_2|OR2|Y~0_combout\,
	datac => \add_instance|fourbit_adder_1|Full_adder_3|XOR2|Y~0_combout\,
	datad => \input_vector[3]~input_o\,
	combout => \add_instance|fourbit_adder_1|Full_adder_3|XOR2|Y~1_combout\);

-- Location: UNVM_X0_Y22_N40
\~QUARTUS_CREATED_UNVM~\ : fiftyfivenm_unvm
-- pragma translate_off
GENERIC MAP (
	addr_range1_end_addr => -1,
	addr_range1_offset => -1,
	addr_range2_offset => -1,
	is_compressed_image => "false",
	is_dual_boot => "false",
	is_eram_skip => "false",
	max_ufm_valid_addr => -1,
	max_valid_addr => -1,
	min_ufm_valid_addr => -1,
	min_valid_addr => -1,
	part_name => "quartus_created_unvm",
	reserve_block => "true")
-- pragma translate_on
PORT MAP (
	nosc_ena => \~QUARTUS_CREATED_GND~I_combout\,
	xe_ye => \~QUARTUS_CREATED_GND~I_combout\,
	se => \~QUARTUS_CREATED_GND~I_combout\,
	busy => \~QUARTUS_CREATED_UNVM~~busy\);

-- Location: ADCBLOCK_X25_Y34_N0
\~QUARTUS_CREATED_ADC1~\ : fiftyfivenm_adcblock
-- pragma translate_off
GENERIC MAP (
	analog_input_pin_mask => 0,
	clkdiv => 1,
	device_partname_fivechar_prefix => "none",
	is_this_first_or_second_adc => 1,
	prescalar => 0,
	pwd => 1,
	refsel => 0,
	reserve_block => "true",
	testbits => 66,
	tsclkdiv => 1,
	tsclksel => 0)
-- pragma translate_on
PORT MAP (
	soc => \~QUARTUS_CREATED_GND~I_combout\,
	usr_pwd => VCC,
	tsen => \~QUARTUS_CREATED_GND~I_combout\,
	chsel => \~QUARTUS_CREATED_ADC1~_CHSEL_bus\,
	eoc => \~QUARTUS_CREATED_ADC1~~eoc\);

-- Location: ADCBLOCK_X25_Y33_N0
\~QUARTUS_CREATED_ADC2~\ : fiftyfivenm_adcblock
-- pragma translate_off
GENERIC MAP (
	analog_input_pin_mask => 0,
	clkdiv => 1,
	device_partname_fivechar_prefix => "none",
	is_this_first_or_second_adc => 2,
	prescalar => 0,
	pwd => 1,
	refsel => 0,
	reserve_block => "true",
	testbits => 66,
	tsclkdiv => 1,
	tsclksel => 0)
-- pragma translate_on
PORT MAP (
	soc => \~QUARTUS_CREATED_GND~I_combout\,
	usr_pwd => VCC,
	tsen => \~QUARTUS_CREATED_GND~I_combout\,
	chsel => \~QUARTUS_CREATED_ADC2~_CHSEL_bus\,
	eoc => \~QUARTUS_CREATED_ADC2~~eoc\);

ww_output_vector(0) <= \output_vector[0]~output_o\;

ww_output_vector(1) <= \output_vector[1]~output_o\;

ww_output_vector(2) <= \output_vector[2]~output_o\;

ww_output_vector(3) <= \output_vector[3]~output_o\;

ww_output_vector(4) <= \output_vector[4]~output_o\;
END structure;


