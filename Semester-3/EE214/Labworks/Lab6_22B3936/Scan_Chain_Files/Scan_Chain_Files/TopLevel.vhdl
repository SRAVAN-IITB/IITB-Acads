library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
				
entity TopLevel is
  port (
    state_out : out std_logic_vector(4 downto 0)
	 --connect to Krypton pins - 58,57,55,53,52
	 --OR MAX X pins - ...
  );
end TopLevel; 

architecture Struct of TopLevel is


  ----------------------------------------------------------------
  --  edit the following lines to set the number of i/o's of your
  --  DUT.
  ----------------------------------------------------------------
  constant number_of_inputs  : integer := 3;  -- # input bits to your design.
  constant number_of_outputs : integer := 2;  -- # output bits from your design.
  ----------------------------------------------------------------
  ----------------------------------------------------------------

   -- declare DUT component
  component DUT  is
     port(input_vector: in std_logic_vector(number_of_inputs-1 downto 0);
        output_vector: out std_logic_vector(number_of_outputs-1 downto 0));
  end component DUT;
  
	-- cross check with jtag_<version>/jtag_inst.vhd
  component v_jtag is
		port (
			virtual_jtag_tdi                : out std_logic;                                       -- tdi
			virtual_jtag_tdo                : in  std_logic                    := 'X';             -- tdo
			virtual_jtag_ir_in              : out std_logic_vector(0 downto 0);                    -- ir_in
			virtual_jtag_ir_out             : in  std_logic_vector(0 downto 0) := (others => 'X'); -- ir_out
			virtual_jtag_virtual_state_cdr  : out std_logic;                                       -- virtual_state_cdr
			virtual_jtag_virtual_state_sdr  : out std_logic;                                       -- virtual_state_sdr
			virtual_jtag_virtual_state_e1dr : out std_logic;                                       -- virtual_state_e1dr
			virtual_jtag_virtual_state_pdr  : out std_logic;                                       -- virtual_state_pdr
			virtual_jtag_virtual_state_e2dr : out std_logic;                                       -- virtual_state_e2dr
			virtual_jtag_virtual_state_udr  : out std_logic;                                       -- virtual_state_udr
			virtual_jtag_virtual_state_cir  : out std_logic;                                       -- virtual_state_cir
			virtual_jtag_virtual_state_uir  : out std_logic;                                       -- virtual_state_uir
			virtual_jtag_tms                : out std_logic;                                       -- tms
			virtual_jtag_jtag_state_tlr     : out std_logic;                                       -- jtag_state_tlr
			virtual_jtag_jtag_state_rti     : out std_logic;                                       -- jtag_state_rti
			virtual_jtag_jtag_state_sdrs    : out std_logic;                                       -- jtag_state_sdrs
			virtual_jtag_jtag_state_cdr     : out std_logic;                                       -- jtag_state_cdr
			virtual_jtag_jtag_state_sdr     : out std_logic;                                       -- jtag_state_sdr
			virtual_jtag_jtag_state_e1dr    : out std_logic;                                       -- jtag_state_e1dr
			virtual_jtag_jtag_state_pdr     : out std_logic;                                       -- jtag_state_pdr
			virtual_jtag_jtag_state_e2dr    : out std_logic;                                       -- jtag_state_e2dr
			virtual_jtag_jtag_state_udr     : out std_logic;                                       -- jtag_state_udr
			virtual_jtag_jtag_state_sirs    : out std_logic;                                       -- jtag_state_sirs
			virtual_jtag_jtag_state_cir     : out std_logic;                                       -- jtag_state_cir
			virtual_jtag_jtag_state_sir     : out std_logic;                                       -- jtag_state_sir
			virtual_jtag_jtag_state_e1ir    : out std_logic;                                       -- jtag_state_e1ir
			virtual_jtag_jtag_state_pir     : out std_logic;                                       -- jtag_state_pir
			virtual_jtag_jtag_state_e2ir    : out std_logic;                                       -- jtag_state_e2ir
			virtual_jtag_jtag_state_uir     : out std_logic;                                       -- jtag_state_uir
			tck_clk                         : out std_logic                                        -- clk
		);
	end component v_jtag;
   
   -- declare I/O signals to DUT component
	signal reg: std_logic_vector(7 downto 0);
	signal data,dut_input: std_logic_vector(number_of_inputs-1 downto 0);
	signal dut_output, datau:  std_logic_vector(number_of_outputs-1 downto 0);
	signal cdr, sdr, udr, e1dr, e2dr, tck ,step1, step2, tms, jsdr, jsir: std_logic; 
	signal io: std_logic:='0';
	
	
	begin
	
	--alternate input and output every cycle
	io_flip: process(tck)
	begin
	if(rising_edge(tck) and udr='1') then
		io <= not io;---1 =output, 0=input
	end if;
	end process;
	
		
	--read data from or write data to internal scan chain depending upon io 
	shift_reg: process(tck)
	begin
	shift_input: if(rising_edge(tck) and sdr='1' and udr='0' and io='0') then
		if(number_of_inputs>1) then
			data<=reg(0) & data(number_of_inputs-1 downto 1);			
		else
			data(0)<=reg(0);
		end if;
	end if;
	
	shift_output: if(rising_edge(tck) and sdr='1' and udr='0' and io='1') then
		if(number_of_outputs>1) then
			datau<="0" & datau(number_of_outputs-1 downto 1) ;	
		else
			datau<="0";
		end if;
	end if;
	
	capture_input: if(rising_edge(tck) and udr='1' and sdr='0' and io ='0') then
		dut_input<= data ;
	end if;
	
	capture_output: if(rising_edge(tck) and cdr='1' and sdr='0' and io ='1') then
		datau<= dut_output ;
	end if;
	
	end process;
		
	--sld instance required for SCAN CHAIN internal logic update the signal names according to Quartus version avoid positional assignment
	u0 : v_jtag
	port map (
		virtual_jtag_tdi                => reg(0),                -- jtag.tdi
		virtual_jtag_tdo                => datau(0),                --     .tdo
		virtual_jtag_ir_in              => reg(7 downto 7),              --     .ir_in
		virtual_jtag_ir_out             => open,             --     .ir_out
		virtual_jtag_virtual_state_cdr  => cdr,  --     .virtual_state_cdr
		virtual_jtag_virtual_state_sdr  => sdr,  --     .virtual_state_sdr
		virtual_jtag_virtual_state_e1dr => e1dr, --     .virtual_state_e1dr
		virtual_jtag_virtual_state_pdr  => open,  --     .virtual_state_pdr
		virtual_jtag_virtual_state_e2dr => e2dr, --     .virtual_state_e2dr
		virtual_jtag_virtual_state_udr  => udr,  --     .virtual_state_udr
		virtual_jtag_virtual_state_cir  => open,  --     .virtual_state_cir
		virtual_jtag_virtual_state_uir  => open,  --     .virtual_state_uir
		virtual_jtag_tms                => tms,                 --  tck.clk
		virtual_jtag_jtag_state_tlr     => open,                                       --     .tms
		virtual_jtag_jtag_state_rti     => open,                                       --     .jtag_state_tlr
		virtual_jtag_jtag_state_sdrs    => open,                                       --     .jtag_state_rti
		virtual_jtag_jtag_state_cdr     => open,                                       --     .jtag_state_sdrs
		virtual_jtag_jtag_state_sdr     => open,                                       --     .jtag_state_cdr
		virtual_jtag_jtag_state_e1dr    => jsdr,                                       --     .jtag_state_sdr
		virtual_jtag_jtag_state_pdr     => open,                                       --     .jtag_state_e1dr
		virtual_jtag_jtag_state_e2dr    => open,                                       --     .jtag_state_pdr
		virtual_jtag_jtag_state_udr     => open,                                       --     .jtag_state_e2dr
		virtual_jtag_jtag_state_sirs    => open,                                       --     .jtag_state_udr
		virtual_jtag_jtag_state_cir     => open,                                       --     .jtag_state_sirs
		virtual_jtag_jtag_state_sir     => open,                                       --     .jtag_state_cir
		virtual_jtag_jtag_state_e1ir    => jsir,                                       --     .jtag_state_sir
		virtual_jtag_jtag_state_pir     => open,                                       --     .jtag_state_e1ir
		virtual_jtag_jtag_state_e2ir    => open,                                       --     .jtag_state_pir
		virtual_jtag_jtag_state_uir     => open,                                       --     .jtag_state_e2ir
		tck_clk                         => tck                                       --     .jtag_state_uir
	);
	
	--student's DUT
	-- connections between DUT and Scan_Chain 
	--use negative edge of UDR as clk is auto clocking is needed
	dut_instance: DUT port map(input_vector=>dut_input, output_vector=>dut_output);

	--just for confirming the correct design is in Krypton
	state_out <= jsdr & cdr & sdr & udr & dut_output(number_of_outputs-1);

end Struct;
