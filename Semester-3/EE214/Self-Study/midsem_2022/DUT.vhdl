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
	Component compare_fourbit is
		port(X, Y : in std_logic_vector(3 downto 0);
			  Z, N, V, L, G : out std_logic);
	end component compare_fourbit;
begin

   -- input/output vector element ordering is critical,
   -- and must match the ordering in the trace file!
   add_instance: compare_fourbit port map (X => input_vector(7 downto 4), Y => input_vector(3 downto 0),
											Z => output_vector(4), N => output_vector(3), V => output_vector(2), 
											L => output_vector(1), G => output_vector(0));

end DutWrap;

