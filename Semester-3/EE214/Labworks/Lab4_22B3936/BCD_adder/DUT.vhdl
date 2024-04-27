-- A DUT entity is used to wrap your design.
--  This example shows how you can do this for the
--  Full-adder.

library ieee;
use ieee.std_logic_1164.all;

entity DUT is
   port(input_vector: in std_logic_vector(7 downto 0);
       	output_vector: out std_logic_vector(4 downto 0));
end entity;

architecture DutWrap of DUT is
	-- Instantiate your own top Module component in place of ALU_1
	Component BCD_Adder is
port (A_3, A_2, A_1, A_0, B_3, B_2, B_1, B_0: in std_logic; Coutt, S_3, S_2, S_1, S_0 : out std_logic);
end  component;
signal A,B: std_logic_vector(3 downto 0);
begin
--	A <= A_3 & A_2 & A_1 & A_0;
--	B <= B_3 & B_2 & B_1 & B_0;
   -- input/output vector element ordering is critical,
   -- and must match the ordering in the trace file!
   add_instance: BCD_Adder port map (A_3 => input_vector(7), A_2 => input_vector(6),
												A_1 => input_vector(5), A_0 => input_vector(4),  
												B_3 => input_vector(3), B_2 => input_vector(2),
												B_1 => input_vector(1), B_0 => input_vector(0),
											Coutt => output_vector(4), S_3 => output_vector(3),
											S_2 => output_vector(2), S_1 => output_vector(1), 
											S_0 => output_vector(0));

end DutWrap;

