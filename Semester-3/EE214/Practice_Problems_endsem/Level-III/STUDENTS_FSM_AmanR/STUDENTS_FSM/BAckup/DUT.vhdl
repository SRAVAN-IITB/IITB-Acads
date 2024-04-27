-- A DUT entity is used to wrap your design.
--  This example shows how you can do this for the
--  Full-adder.

library ieee;
use ieee.std_logic_1164.all;

entity DUT is
   port(input_vector: in std_logic_vector(6 downto 0);
       	output_vector: out std_logic_vector(0 downto 0));
end entity;

architecture DutWrap of DUT is
	-- Instantiate your own top Module component in place of ALU_1
	component STUDENTS_FSM is
		port ( A,B,C,D,E,reset, clock : in std_logic; output : out std_logic);
	end component STUDENTS_FSM;
begin

   -- input/output vector element ordering is critical,
   -- and must match the ordering in the trace file!
   add_instance: STUDENTS_FSM port map (A => input_vector(6),
												B => input_vector(5),
												C => input_vector(4),
												D => input_vector(3),
												E => input_vector(2),
												reset => input_vector(1),
												clock => input_vector(0),
											output => output_vector(0));

end DutWrap;

