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
	Component Fourbit_Adder is
		generic (operand_width : integer := 4);
		port(A, B : in std_logic_vector(operand_width-1 downto 0);
			  Op   : out std_logic_vector(4 downto 0));
	end component Fourbit_Adder;
begin

   -- input/output vector element ordering is critical,
   -- and must match the ordering in the trace file!
   add_instance: Fourbit_Adder port map (A => input_vector(7 downto 4), B => input_vector(3 downto 0),
											Op => output_vector);

end DutWrap;

