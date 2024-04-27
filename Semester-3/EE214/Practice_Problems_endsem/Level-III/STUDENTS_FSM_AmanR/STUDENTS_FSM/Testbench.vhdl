library ieee;
use ieee.std_logic_1164.all;

entity Testbench is
end entity;

architecture TB_BHV of Testbench is

	component STUDENTS_FSM is
		port ( A,B,C,D,E,reset, clock : in std_logic; output : out std_logic);
	end component STUDENTS_FSM;
	
	signal A,B,C,D,E, reset, clock : std_logic := '0';
	signal output : std_logic;
	
begin

	clock <= not clock after 10ns;
	A <= not A after 9ns;
	B <= not B after 8ns;
	C <= not C after 7ns;
	D <= not D after 6ns;
	E <= not E after 4ns;
	
	FSM1: STUDENTS_FSM port map( A,B,C,D,E,reset, clock, output);
	
end_sim: process
begin
	wait for 5000ns;
	assert false
	report "End Simulation"
	severity failure;
end process end_sim;
	
end TB_BHV;