-- A DUT entity is used to wrap your design so that we can combine it with testbench.

library IEEE;
use IEEE.std_logic_1164.all;

entity DUT is
    port(input_vector: in std_logic_vector(17 downto 0);
       	output_vector: out std_logic_vector(15 downto 0));
end entity;

architecture DutWrap of DUT is
   component Reg_16BIT is
		port (Write_EN, clk	: 	in std_logic;
				data_in			:	in std_logic_vector(15 downto 0);
				data_out			: out std_logic_vector(15 downto 0)
				);
	end component;

begin
   add_instance: Reg_16BIT 
			port map(Write_EN => input_vector(17), 
						clk 		=> input_vector(16),
						data_in  =>	input_vector(15 downto 0),
						data_out	=> output_vector(15 downto 0)
				);
end DutWrap;