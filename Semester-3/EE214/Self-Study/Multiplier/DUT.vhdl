-- A DUT entity is used to wrap your design.
--  This example shows how you can do this for the
--  Full-adder.

library ieee;
use ieee.std_logic_1164.all;

entity DUT is
   port(input_vector: in std_logic_vector(6 downto 0);
        output_vector: out std_logic_vector(6 downto 0));
end entity;

architecture DutWrap of DUT is
	-- Instantiate your own top Module component
	Component mult is
		port (a: in std_logic_vector( 3 downto 0);
		b: in std_logic_vector( 2 downto 0); 
		o: out std_logic_vector(6 downto 0));
	end component mult;
begin

   -- input/output vector element ordering is critical,
   -- and must match the ordering in the trace file!
   add_instance: mult port map (a => input_vector(6 downto 3), b => input_vector(2 downto 0),
										  o => output_vector(6 downto 0));

end DutWrap;

