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
			Q, NotQ : out std_logic);
end component JK_FlipFlop;

signal J, K, Clock, Enable, Preset, Reset, Q, NotQ : std_logic;

begin
	FlipFlop : JK_FlipFlop port map (J => J, K => K,
									  Clock => Clock, 
									  Enable => Enable,
									  Preset => Preset, Reset => Reset,
									  Q => Q, NotQ => NotQ);
									  
	ClockProcess : process
	begin
		Clock <= '1';
		wait for 10 ns;
		Clock <= '0';
		wait for 10 ns;
	end process;
	
	Preset <= '1', '0' after 10 ns;
	Reset <= '1', '0' after 20 ns;
	Enable <= '0', '1' after 30 ns;
	
	J <= '1', '0' after 40 ns, '0' after 50 ns, '1' after 60 ns;
	K <= '0', '1' after 40 ns, '0' after 50 ns, '1' after 60 ns;
	
	end_sim : process
	begin
		wait for 150ns;
		ASSERT false
		REPORT "Sim Ended"
		severity failure;
	end process end_sim;
		
end TB_BHV;
