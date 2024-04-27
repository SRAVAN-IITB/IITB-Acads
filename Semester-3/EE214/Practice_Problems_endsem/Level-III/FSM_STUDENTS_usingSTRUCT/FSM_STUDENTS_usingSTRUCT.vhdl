library IEEE;
use IEEE.std_logic_1164.all;
use IEEE.std_logic_arith.all;

entity FSM_STUDENTS_usingSTRUCT is
	port( inp : in std_logic_vector(4 downto 0);
			reset, clock : in std_logic;
			outp : out std_logic);
end FSM_STUDENTS_usingSTRUCT;

architecture STRUCT of FSM_STUDENTS_usingSTRUCT is
begin
	---------------Define state type here-----------------------------
	type state is (rst, s1, s2, s3, s4, s5, s6, s7);
		
	
	
end STRUCT;