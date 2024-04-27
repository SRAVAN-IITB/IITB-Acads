library IEEE;
use IEEE.std_logic_1164.all;

entity Testbench is
end entity;

architecture TB_BHV of Testbench is
	component Synch_Down_Ctr_4bit is
		PORT(clock : in std_logic;
				counter : out std_logic_vector(3 DOWNTO 0));
	end component Synch_Down_Ctr_4bit;
	
	signal clock : std_logic := '0';
	signal counter : std_logic_vector(3 DOWNTO 0);
	
	begin
		clock <= not clock after 10ns;
		
		Synch_Down_Ctr_4bit_1 : Synch_Down_Ctr_4bit PORT MAP(clock, counter);
		
		end_sim : PROCESS
			begin
				wait for 150ns;
				ASSERT False
				REPORT "Sim ended"
				SEVERITY failure;
		end PROCESS end_sim;
				
end TB_BHV;