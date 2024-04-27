library IEEE;
use IEEE.std_logic_1164.all;

entity Testbench is
end entity;

architecture TB_BHV of Testbench is
component JK_FlipFlop is
	port (J, K : in std_logic;
			Clock : in std_logic;
			Enable : in std_logic;
			Preset, Reset : in std_logic;
			Q, QBar : out std_logic);
end component JK_FlipFlop;

signal J, K, Clock, Enable, Preset, Reset: std_logic := '0';
signal Q, QBar : std_logic;

begin
	FlipFlop : JK_FlipFlop port map (J => J, K => K,
									  Clock => Clock, 
									  Enable => Enable,
									  Preset => Preset, Reset => Reset,
									  Q => Q, QBar => QBar);
	
	clock <= not clock after 10ns;
	Preset <= '1', '0' after 10 ns;
	Reset <= '1', '0' after 20 ns;
	Enable <= '0', '1' after 30 ns;
	J <= '0', '0' after 20 ns, '1' after 30 ns, '1' after 40 ns;
	K <= '0', '1' after 20 ns, '0' after 30 ns, '1' after 40 ns;
	
	end_sim : process
	begin
		wait for 400ns;
		ASSERT False
		REPORT "Sim ended"
		SEVERITY failure;
	end process end_sim;
		
end TB_BHV;
