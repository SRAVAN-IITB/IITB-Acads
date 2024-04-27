-- A DUT entity is used to wrap your design so that we can combine it with testbench.
-- This example shows how you can do this for the OR Gate

library ieee;
use ieee.std_logic_1164.all;

entity DUT is
    port(input_vector: in std_logic_vector(2 downto 0);
       	output_vector: out std_logic_vector(0 downto 0));
end entity;

architecture DutWrap of DUT is
   component MUX_2x1 is
     port(in0,in1,S: in std_logic; Y0: out std_logic);
   end component;
begin

   -- input/output vector element ordering is critical,
   -- and must match the ordering in the trace file!
   add_instance: MUX_2x1 
			port map (
					-- order of inputs B A
					in1 => input_vector(2),
					in0=> input_vector(1),
					S=> input_vector(0),
               -- order of output OUTPUT
					Y0 => output_vector(0));
end DutWrap;