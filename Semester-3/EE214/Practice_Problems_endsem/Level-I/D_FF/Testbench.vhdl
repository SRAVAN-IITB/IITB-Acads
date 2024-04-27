library IEEE;
use IEEE.std_logic_1164.all;

entity Testbench is
end entity;

architecture TB_BHV of Testbench is

	component D_FF is
		port(clock, D, DBar, Reset, Preset, Enable : in std_logic; Q, QBar : out std_logic);
	end component D_FF; 
	signal clock, D, DBar, Reset, Preset, Enable : std_logic := '0'; 
	signal Q, QBar : std_logic;
	
	begin
		clock <= not clock after 10ns;
		D <= not D after 8ns;
		DBar <= not DBar after 6ns;
		Preset <= '1', '0' after 15ns;
		Reset <= '1', '0' after 30ns;
		Enable <= '1';
	D_FF1: D_FF PORT MAP(clock, D, Dbar, Reset, Preset, Enable, Q, QBar);
	
	end_sim : process
	begin
		wait for 150ns;
		ASSERT false
		REPORT "Sim Ended"
		severity failure;
	end process end_sim;
	
end TB_BHV;
