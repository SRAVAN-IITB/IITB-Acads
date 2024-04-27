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
component OR_gate is
   port (A, B: in std_logic; Y: out std_logic);
  end component OR_gate;

begin

   -- input/output vector element ordering is critical,
   -- and must match the ordering in the trace file!
   add_instance: OR_gate port map (A => input_vector(1), B => input_vector(0),
											Y => output_vector(0));
end DutWrap;