library IEEE;
use IEEE.std_logic_1164.all;
use IEEE.numeric_std.all;

entity Testbench is
end entity;

architecture Struct of Testbench is
    component Register6bit is
        port(Clk        : in  std_logic;
			 Reset      : in  std_logic;
             data_in    : in  std_logic_vector(5 downto 0);
             data_out   : out std_logic_vector(5 downto 0));
    end component;

    signal Clk, Reset : std_logic;
    signal data_in    : std_logic_vector(5 downto 0);
    signal data_out   : std_logic_vector(5 downto 0);
    signal i, j       : integer;

begin
    DUT: Register6bit port map (Clk => Clk, Reset => Reset, data_in => data_in, data_out => data_out);

    Reset <= '1', '0' after 40ns;

    Prc1 : process
    begin
        for j in 0 to 63 loop
            Clk <= '0';
            wait for 10ns;
            Clk <= '1';
            wait for 10ns;
        end loop;
    end process Prc1 ;

    Prc2: process
    begin
        for i in 0 to 63 loop
            data_in <= std_logic_vector(to_unsigned(i,6));
            wait for 20ns;
        end loop;
    end process Prc2;
end architecture Struct;