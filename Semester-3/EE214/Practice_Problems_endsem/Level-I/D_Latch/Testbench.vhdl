library IEEE;
use IEEE.std_logic_1164.all;

entity Testbench is
end Testbench;

architecture TB_BHV of Testbench is
	component D_Latch is 
		port (clock, D :in std_logic ; Q, QBAR : out std_logic);
	end component D_Latch;
	signal clock, D : std_logic := '0';
	signal Q, QBAR : std_logic;
begin
	clock <= not clock after 10ns;
	D <= not D after 15ns;
	
	D_L1: D_Latch port map(clock, D, Q, QBAR);
	
	end_sim: process
	begin
		wait for 150ns;
		assert false
		report "Sim Ended"
		severity failure;
	end process end_sim;
	
end TB_BHV;