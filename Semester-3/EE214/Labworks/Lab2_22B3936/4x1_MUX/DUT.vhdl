-- A DUT entity is used to wrap your design so that we can combine it with testbench.
-- This example shows how you can do this for the OR Gate

library ieee;
use ieee.std_logic_1164.all;

entity DUT is
    port(input_vector: in std_logic_vector(5 downto 0);
       	output_vector: out std_logic_vector(0 downto 0));
end entity;

architecture DutWrap of DUT is
   component MUX_4x1 is
		port ( In4, In3, In2, In1, S2, S1: in std_logic; Y0: out std_logic);
   end component;
begin

   -- input/output vector element ordering is critical,
   -- and must match the ordering in the trace file!
   add_instance: MUX_4x1 
			port map (
					-- order of inputs B A
					In4 => input_vector(5),
					In3 => input_vector(4),
					In2 => input_vector(3),
					In1 => input_vector(2),
					S2=> input_vector(1),
					S1=> input_vector(0),
               -- order of output OUTPUT
					Y0 => output_vector(0));
end DutWrap;