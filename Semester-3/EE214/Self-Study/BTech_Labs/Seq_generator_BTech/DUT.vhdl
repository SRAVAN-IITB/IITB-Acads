-- A DUT entity is used to wrap your design.
--  This example shows how you can do this for the
--  Full-adder.

library ieee;
use ieee.std_logic_1164.all;

entity DUT is
   port(input_vector: in std_logic_vector(1 downto 0);
       	output_vector: out std_logic_vector(0 downto 0));
end entity;

architecture DutWrap of DUT is
	-- Instantiate your own top Module component in place of Sequence_generator_stru_dataflow
	Component Sequence_generator_stru_dataflow is
		port (reset,clock: in std_logic;
		y:out std_logic);
	end component Sequence_generator_stru_dataflow;
begin

   -- input/output vector element ordering is critical,
   -- and must match the ordering in the trace file!
   add_instance: Sequence_generator_stru_dataflow port map (
					-- order of inputs B A C
					reset => input_vector(1),
					clock => input_vector(0),
               -- order of output OUTPUT
					y => output_vector(0));

end DutWrap;

