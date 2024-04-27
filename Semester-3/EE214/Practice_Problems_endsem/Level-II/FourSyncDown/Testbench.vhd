library ieee;
use ieee.std_logic_1164.all;

entity Testbench is
end Testbench;

architecture TB_BHV of Testbench is
	component FourSyncDown is
	port (Clock, Enable, Preset, Reset : in std_logic;
			Counter : out std_logic_vector(3 downto 0));
	end component FourSyncDown;
	
	signal Clock : std_logic := '0';
	signal Enable, Preset, Reset : std_logic;
	signal Counter : std_logic_vector(3 downto 0);
	
	begin
--	counter <= "1111";
--	wait for 10ns;
	CounterInst : FourSyncDown port map (Clock => Clock, Enable => Enable,
													 Preset => Preset, Reset => Reset,
													 Counter => Counter);
	 
	Preset <= '1', '0' after 30ns;
	Reset <= '1', '0' after 15ns;
	clock <= not clock after 10ns;
--		D <= not D after 8ns;
		Enable <= '1';
	
	end_sim : process
	begin
		wait for 400ns;
		ASSERT False
		REPORT "Sim ended"
		SEVERITY failure;
	end process end_sim;
end TB_BHV;