LIBRARY ieee;
USE ieee.std_logic_1164.ALL;
 
ENTITY tb_bcd_7seg IS
END tb_bcd_7seg;
 
ARCHITECTURE behavior OF tb_bcd_7seg IS
 
-- Component Declaration for the Unit Under Test (UUT)
 
COMPONENT BCDtoSevenSegment_Display_bhv_usingCASE
	PORT(
	BCDin : IN std_logic_vector(3 downto 0);
	Seven_Segment : OUT std_logic_vector(6 downto 0)
	);
END COMPONENT;
 
--Inputs
signal BCDin : std_logic_vector(3 downto 0) := (others => '0');
 
--Outputs
signal Seven_Segment : std_logic_vector(6 downto 0);
 
BEGIN
-- Instantiate the Unit Under Test (UUT)
	uut: BCDtoSevenSegment_Display_bhv_usingCASE PORT MAP (
		BCDin => BCDin,
		Seven_Segment => Seven_Segment
	);
 
	-- Stimulus process
	stim_proc: process
	begin
	 
		BCDin <= "0000";
		wait for 100 ns;
		
		BCDin <= "0001";
		wait for 100 ns;
		
		BCDin <= "0010";
		wait for 100 ns;
		
		BCDin <= "0011";
		wait for 100 ns;
		
		BCDin <= "0100";
		wait for 100 ns;
		
		BCDin <= "0101";
		wait for 100 ns;
		
		BCDin <= "0110";
		wait for 100 ns;
		
		BCDin <= "0111";
		wait for 100 ns;
		
		BCDin <= "1000";
		wait for 100 ns;
		
		BCDin <= "1001";
		wait for 100 ns;
		
	end process;
 
END behavior;
