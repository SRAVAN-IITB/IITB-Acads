-- A DUT entity is used to wrap your design so that we can combine it with testbench.
-- This example shows how you can do this for the OR Gate

library ieee;
use ieee.std_logic_1164.all;

entity DUT is
    port(input_vector: in std_logic_vector(1 downto 0);
       	output_vector: out std_logic_vector(0 downto 0));
end entity;

architecture DutWrap of DUT is
   component SequenceGenerator is
		port (reset,clock: in std_logic; y:out std_logic);
	end component SequenceGenerator;
begin

   -- input/output vector element ordering is critical,
   -- and must match the ordering in the trace file!
   add_instance: SequenceGenerator 
			port map (
					-- order of inputs B A
					reset => input_vector(1),
					clock => input_vector(0),
               -- order of output OUTPUT
					y => output_vector(0));
end DutWrap;