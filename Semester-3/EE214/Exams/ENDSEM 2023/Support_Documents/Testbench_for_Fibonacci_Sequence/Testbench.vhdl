library IEEE;
use IEEE.std_logic_1164.all;
use IEEE.numeric_std.all;

entity Testbench is
end entity;

architecture Struct of Testbench is
    component FibGenerator is
       port(Clk, Reset : in std_logic;
            Seq_out    : out std_logic_vector(5 downto 0));
    end component;

    component Clock_Divider is
        port ( Clk_50  : in std_logic;
               Reset   : in std_logic;
               Clk_Out : out std_logic );
    end component;

    signal Clk, Reset : std_logic;
    signal Clk_Out    : std_logic;
    signal Seq_out    : std_logic_vector(5 downto 0);
    signal i, j       : integer;
    --signal a, b, c    : std_logic_vector(5 downto 0);
    signal golden_out : std_logic_vector(5 downto 0);
    signal flag       : std_logic;
begin
    DUT: FibGenerator port map (Clk => Clk, Reset => Reset, Seq_out => Seq_out);
    CLK1 : Clock_Divider port map (Clk, Reset, Clk_Out);

    Reset <= '1', '0' after 1ms;

    Prc1 : process
    begin
        Clk <= '0';
        wait for 10ns;
        Clk <= '1';
        wait for 10ns;
    end process Prc1 ;

    Prc2: process (Clk_Out)
        variable a    : std_logic_vector(6 downto 0) := "0000000"; 
        variable b, c : std_logic_vector(6 downto 0) := "0000001";
    begin
        if (rising_edge(Clk_Out)) then
            c := a;
            a := std_logic_vector(unsigned(a) + unsigned(b));
            b := c;
            if (a(6) = '1') then
                a := "0000001";
                b := "0000000";
            end if;
            golden_out <= a(5 downto 0);
        end if;
    end process Prc2;

    flag <= '1' when golden_out = Seq_out else '0';
end architecture Struct;
