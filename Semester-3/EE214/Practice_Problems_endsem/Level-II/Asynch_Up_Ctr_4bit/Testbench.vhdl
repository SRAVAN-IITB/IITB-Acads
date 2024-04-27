library IEEE;
use IEEE.std_logic_1164.all;

entity Testbench is
end entity Testbench;

architecture TB_BHV of Testbench is

	component Asynch_Up_Ctr_4bit is
	port (Clock, Reset : in std_logic;
			Counter : out std_logic_vector(3 downto 0));
	end component Asynch_Up_Ctr_4bit;
	
	signal Reset : std_logic;
	signal clock : std_logic := '0';
--	signal Q : std_logic_vector(3 downto 0) ;
--	signal QBar : std_logic_vector(3 downto 0) ;	
	signal Counter : std_logic_vector(3 downto 0);
	
	begin
			clock <= not clock after 10ns;

			Reset <= '1', '0' after 5 ns;
--			
		Asynch_Up_Ctr_4bit1 : Asynch_Up_Ctr_4bit PORT MAP(clock, Reset, Counter);
--		Asynch_Up_Ctr_4bit1 : Asynch_Up_Ctr_4bit PORT MAP(Clock => Clock, Counter => Counter, Enable => Enable, 
--																			Preset => Preset, Reset => Reset);
		
		end_sim : process
		begin
			wait for 400ns;
			ASSERT False
			REPORT "Sim ended"
			SEVERITY failure;
		end process end_sim;
		
end TB_BHV;