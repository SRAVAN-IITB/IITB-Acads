library ieee;
use ieee.std_logic_1164.all;

entity Testbench is
end entity Testbench;

architecture TB_BHV of Testbench is
	component FourBitAsyncUp is
	port (Clock : in std_logic;
			Counter : out std_logic_vector(3 downto 0);
			Enable, Preset, Reset : in std_logic);
	end component FourBitAsyncUp;
	
	signal Clock, Enable, Preset, Reset : std_logic := '0';
	signal Counter : std_logic_vector(3 downto 0);
	
	begin
	CounterInst : FourBitAsyncUp port map (Clock => Clock, Counter => Counter, Enable => Enable, 
															Preset => Preset, Reset => Reset);
		clock <= not clock after 10ns;
		Preset <= '1', '0' after 15ns;
		Reset <= '1', '0' after 30ns;
		Enable <= '1';
	
	end_sim : process
	begin
		wait for 400ns;
		ASSERT False
		REPORT "Sim ended"
		SEVERITY failure;
	end process end_sim;
		
end TB_BHV;