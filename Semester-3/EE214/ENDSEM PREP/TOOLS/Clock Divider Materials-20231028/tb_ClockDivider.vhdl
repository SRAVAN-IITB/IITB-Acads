library ieee;
use ieee.std_logic_1164.all;

entity tb_ClockDivider is
end tb_ClockDivider;

architecture ClockTB of tb_ClockDivider is
	
	-- Component Instaiation of the Entity you want to verify using testbench
	component Clock_Divider is
		port ( Clk_50  : in std_logic;
			   Reset   : in std_logic;
			   Clk_Out : out std_logic );
	end component Clock_Divider;
	
	signal Clk_50, Reset, Clk_Out : std_logic;
	
	begin
		-- Port Mapping of the component being instantiated
		ClkDiv1 : Clock_Divider port map (Clk_50 => Clk_50, Reset => Reset, Clk_Out => Clk_Out);
		
		Reset <= '1', '0' after 20 ms;
		
		L1: process  -- In Testbench Process statement does not have sensitivity list
					 -- the Statement written inside process block in testbench will run in a infinite loop
				begin
					Clk_50 <= '0';
					wait for 100 us; -- 100 us is used as Clk_50 freq = 50 MHz, so T = 200 us 
									 -- So T/2 Clk_50 will be OFF and for next T/2 Clk_50 will ON 
					Clk_50 <= '1';
					wait for 100 us;
			end process;
		
end ClockTB;