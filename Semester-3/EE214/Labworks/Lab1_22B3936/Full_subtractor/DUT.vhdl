-- A DUT entity is used to wrap your design.
--  This example shows how you can do this for the
--  Full-adder.

library ieee;
use ieee.std_logic_1164.all;

entity DUT is
   port(input_vector: in std_logic_vector(2 downto 0);
       	output_vector: out std_logic_vector(1 downto 0));
end entity;

architecture DutWrap of DUT is
component Full_subtractor is
   port (A, B, C: in std_logic; Dif, Bor: out std_logic);
  end component Full_subtractor;

begin

   -- input/output vector element ordering is critical,
   -- and must match the ordering in the trace file!
   add_instance: Full_subtractor port map (A => input_vector(2), B => input_vector(1), C => input_vector(0),
											Dif => output_vector(1), Bor => output_vector(0));
end DutWrap;