library IEEE;
use IEEE.std_logic_1164.all;
use ieee.numeric_std.all;

entity Reg_16BIT is
	port (Write_EN, clk	: 	in std_logic;
			data_in			:	in std_logic_vector(15 downto 0);
			data_out			: out std_logic_vector(15 downto 0)
			);
end entity Reg_16BIT;

architecture bhv of Reg_16BIT is
	signal reg_store: std_logic_vector(15 downto 0):="0000000000000000";
	begin
		data_out(15 downto 0) <= reg_store(15 downto 0);
		
		clkProc: process(clk, Write_EN, reg_store)
		begin
			if(clk='1' and clk'event and Write_EN = '1') then
				reg_store(15 downto 0) <= data_in(15 downto 0);
			else
				reg_store(15 downto 0) <= reg_store(15 downto 0);
			end if;
		end process;
		
end architecture bhv;