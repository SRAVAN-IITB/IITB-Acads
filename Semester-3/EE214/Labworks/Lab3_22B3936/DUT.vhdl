-- A DUT entity is used to wrap your design.
--  This example shows how you can do this for the
--  Full-adder.

library ieee;
use ieee.std_logic_1164.all;

entity DUT is
   port(input_vector: in std_logic_vector(4 downto 0);
       	output_vector: out std_logic_vector(0 downto 0));
end entity;

architecture DutWrap of DUT is
	-- Instantiate your own top Module component in place of fibonacci
	Component fibonacci is
		port (in4, in3, in2, in1, in0: in std_logic; Y0: out std_logic);
	end component;
begin

   -- input/output vector element ordering is critical,
   -- and must match the ordering in the trace file!
   add_instance: fibonacci port map (in4 => input_vector(4), in3 => input_vector(3), 
													in2 => input_vector(2), in1 => input_vector(1), 
													in0 => input_vector(0),
											Y0 => output_vector(0));

end DutWrap;

