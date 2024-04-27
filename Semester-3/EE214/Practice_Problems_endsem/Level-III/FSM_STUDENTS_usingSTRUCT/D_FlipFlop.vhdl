library IEEE;
use IEEE.std_logic_1164.all;

entity D_FlipFlop is
	port (D : in std_logic;
			Clock : in std_logic;
			Enable : in std_logic;
			Reset, Preset : in std_logic;
			Q, QBar : out std_logic);
end entity D_FlipFlop;

architecture Struct of D_FlipFlop is

component JK_FlipFlop is
	port (J, K : in std_logic;
			Clock : in std_logic;
			Enable : in std_logic;
			Preset, Reset : in std_logic;
			Q, QBar : out std_logic);
end component JK_FlipFlop;

signal Output, NotOutput, NotD : std_logic;

begin
	NotD <= not D;
	JKFF_Instance : JK_FlipFlop port map (J => D, K => NotD, 
													  Clock => Clock,
													  Enable => Enable,
													  Preset => Preset, 
													  Reset => Reset,
													  Q => Output,
													  QBar => NotOutput);										 
	Q <= Output;
	QBar <= NotOutput;
end Struct;