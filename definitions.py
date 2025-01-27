# important parameters
import math

design_offset = 88 * math.pi / 180 # required to determine design parameters
earth_offset = math.pi # earth offset (180* against the sun)
pi2 = math.pi*2 # circle constant

# Important parameters and statistics for mandala wheel calculation
hex_width = pi2/64
line_width = hex_width/6
color_width = line_width/6
tone_width = color_width/6
base_width = tone_width/5

# Credits to Ra Uru Hu for receiving this info
iching_map = [36, 22, 63, 37, 55, 30, 49, 13, 19, 41, 60, 61, 54, 38, 58, 10, 11, 26, 5, 9, 34, 14, 43, 1, 44, 28, 50, 32, 57, 48, 18, 46, 6, 47, 64, 40, 59, 29, 4, 7, 33, 31, 56, 62, 53, 39, 52, 15, 12, 45, 35, 16, 20, 8, 23, 2, 24, 27, 3, 42, 51, 21, 17, 25]
iching_map.reverse() 

iching = [
    {"hex":1,"hex_font":"䷀","trad_chinese":"乾","pinyin":"qián","english":"Initiating","binary":111111,"od":"02"},
    {"hex":2,"hex_font":"䷁","trad_chinese":"坤","pinyin":"kūn","english":"Responding","binary":"000000","od":"01"},
    {"hex":3,"hex_font":"䷂","trad_chinese":"屯","pinyin":"zhūn","english":"Beginning","binary":"010001","od":50},
    {"hex":4,"hex_font":"䷃","trad_chinese":"蒙","pinyin":"méng","english":"Childhood","binary":100010,"od":49},
    {"hex":5,"hex_font":"䷄","trad_chinese":"需","pinyin":"xū","english":"Needing","binary":"010111","od":35},
    {"hex":6,"hex_font":"䷅","trad_chinese":"訟","pinyin":"sòng","english":"Contention","binary":111010,"od":36},
    {"hex":7,"hex_font":"䷆","trad_chinese":"師","pinyin":"shī","english":"Multitude","binary":"000010","od":13},
    {"hex":8,"hex_font":"䷇","trad_chinese":"比","pinyin":"bǐ","english":"Union","binary":"010000","od":14},
    {"hex":9,"hex_font":"䷈","trad_chinese":"小畜","pinyin":"xiǎochù","english":"LittleAccumulation","binary":110111,"od":16},
    {"hex":10,"hex_font":"䷉","trad_chinese":"履","pinyin":"lǚ","english":"Fulfillment","binary":111011,"od":15},
    {"hex":11,"hex_font":"䷊","trad_chinese":"泰","pinyin":"tài","english":"Advance","binary":"000111","od":12},
    {"hex":12,"hex_font":"䷋","trad_chinese":"否","pinyin":"pǐ","english":"Hindrance","binary":111000,"od":11},
    {"hex":13,"hex_font":"䷌","trad_chinese":"同人","pinyin":"tóngrén","english":"SeekingHarmony","binary":111101,"od":"07"},
    {"hex":14,"hex_font":"䷍","trad_chinese":"大有","pinyin":"dàyǒu","english":"GreatHarvest","binary":101111,"od":"08"},
    {"hex":15,"hex_font":"䷎","trad_chinese":"謙","pinyin":"qiān","english":"Humbleness","binary":"000100","od":10},
    {"hex":16,"hex_font":"䷏","trad_chinese":"豫","pinyin":"yù","english":"Delight","binary":"001000","od":"09"},
    {"hex":17,"hex_font":"䷐","trad_chinese":"隨","pinyin":"suí","english":"Following","binary":"011001","od":18},
    {"hex":18,"hex_font":"䷑","trad_chinese":"蠱","pinyin":"gǔ","english":"Remedying","binary":100110,"od":17},
    {"hex":19,"hex_font":"䷒","trad_chinese":"臨","pinyin":"lín","english":"Approaching","binary":"000011","od":33},
    {"hex":20,"hex_font":"䷓","trad_chinese":"觀","pinyin":"guān","english":"Watching","binary":110000,"od":34},
    {"hex":21,"hex_font":"䷔","trad_chinese":"噬嗑","pinyin":"shìkè","english":"Eradicating","binary":101001,"od":48},
    {"hex":22,"hex_font":"䷕","trad_chinese":"賁","pinyin":"bì","english":"Adorning","binary":100101,"od":47},
    {"hex":23,"hex_font":"䷖","trad_chinese":"剝","pinyin":"bō","english":"FallingAway","binary":100000,"od":43},
    {"hex":24,"hex_font":"䷗","trad_chinese":"復","pinyin":"fù","english":"TurningBack","binary":"000001","od":44},
    {"hex":25,"hex_font":"䷘","trad_chinese":"無妄","pinyin":"wúwàng","english":"WithoutFalsehood","binary":111001,"od":46},
    {"hex":26,"hex_font":"䷙","trad_chinese":"大畜","pinyin":"dàchù","english":"GreatAccumulation","binary":100111,"od":45},
    {"hex":27,"hex_font":"䷚","trad_chinese":"頤","pinyin":"yí","english":"Nourishing","binary":100001,"od":28},
    {"hex":28,"hex_font":"䷛","trad_chinese":"大過","pinyin":"dàguò","english":"GreatExceeding","binary":"011110","od":27},
    {"hex":29,"hex_font":"䷜","trad_chinese":"坎","pinyin":"kǎn","english":"Darkness","binary":"010010","od":30},
    {"hex":30,"hex_font":"䷝","trad_chinese":"離","pinyin":"lí","english":"Brightness","binary":101101,"od":29},
    {"hex":31,"hex_font":"䷞","trad_chinese":"咸","pinyin":"xián","english":"MutualInfluence","binary":"011100","od":41},
    {"hex":32,"hex_font":"䷟","trad_chinese":"恆","pinyin":"héng","english":"LongLasting","binary":"001110","od":42},
    {"hex":33,"hex_font":"䷠","trad_chinese":"遯","pinyin":"dùn","english":"Retreat","binary":111100,"od":19},
    {"hex":34,"hex_font":"䷡","trad_chinese":"大壯","pinyin":"dàzhuàng","english":"GreatStrength","binary":"001111","od":20},
    {"hex":35,"hex_font":"䷢","trad_chinese":"晉","pinyin":"jìn","english":"ProceedingForward","binary":101000,"od":"05"},
    {"hex":36,"hex_font":"䷣","trad_chinese":"明夷","pinyin":"míngyí","english":"BrillianceInjured","binary":"000101","od":"06"},
    {"hex":37,"hex_font":"䷤","trad_chinese":"家人","pinyin":"jiārén","english":"Household","binary":110101,"od":40},
    {"hex":38,"hex_font":"䷥","trad_chinese":"睽","pinyin":"kuí","english":"Diversity","binary":101011,"od":39},
    {"hex":39,"hex_font":"䷦","trad_chinese":"蹇","pinyin":"jiǎn","english":"Hardship","binary":"010100","od":38},
    {"hex":40,"hex_font":"䷧","trad_chinese":"解","pinyin":"xiè","english":"Relief","binary":"001010","od":37},
    {"hex":41,"hex_font":"䷨","trad_chinese":"損","pinyin":"sǔn","english":"Decreasing","binary":100011,"od":31},
    {"hex":42,"hex_font":"䷩","trad_chinese":"益","pinyin":"yì","english":"Increasing","binary":110001,"od":32},
    {"hex":43,"hex_font":"䷪","trad_chinese":"夬","pinyin":"guài","english":"Eliminating","binary":"011111","od":23},
    {"hex":44,"hex_font":"䷫","trad_chinese":"姤","pinyin":"gòu","english":"Encountering","binary":111110,"od":24},
    {"hex":45,"hex_font":"䷬","trad_chinese":"萃","pinyin":"cuì","english":"BringingTogether","binary":"011000","od":26},
    {"hex":46,"hex_font":"䷭","trad_chinese":"升","pinyin":"shēng","english":"GrowingUpward","binary":"000110","od":25},
    {"hex":47,"hex_font":"䷮","trad_chinese":"困","pinyin":"kùn","english":"Exhausting","binary":"011010","od":22},
    {"hex":48,"hex_font":"䷯","trad_chinese":"井","pinyin":"jǐng","english":"Replenishing","binary":"010110","od":21},
    {"hex":49,"hex_font":"䷰","trad_chinese":"革","pinyin":"gé","english":"AbolishingTheOld","binary":"011101","od":"04"},
    {"hex":50,"hex_font":"䷱","trad_chinese":"鼎","pinyin":"dǐng","english":"EstablishingTheNew","binary":101110,"od":"03"},
    {"hex":51,"hex_font":"䷲","trad_chinese":"震","pinyin":"zhèn","english":"TakingAction","binary":"001001","od":57},
    {"hex":52,"hex_font":"䷳","trad_chinese":"艮","pinyin":"gèn","english":"KeepingStill","binary":100100,"od":58},
    {"hex":53,"hex_font":"䷴","trad_chinese":"漸","pinyin":"jiàn","english":"DevelopingGradually","binary":110100,"od":54},
    {"hex":54,"hex_font":"䷵","trad_chinese":"歸妹","pinyin":"guīmèi","english":"MarryingMaiden","binary":"001011","od":53},
    {"hex":55,"hex_font":"䷶","trad_chinese":"豐","pinyin":"fēng","english":"Abundance","binary":"001101","od":59},
    {"hex":56,"hex_font":"䷷","trad_chinese":"旅","pinyin":"lǚ","english":"Travelling","binary":101100,"od":60},
    {"hex":57,"hex_font":"䷸","trad_chinese":"巽","pinyin":"xùn","english":"ProceedingHumbly","binary":110110,"od":51},
    {"hex":58,"hex_font":"䷹","trad_chinese":"兌","pinyin":"duì","english":"Joyful","binary":"011011","od":52},
    {"hex":59,"hex_font":"䷺","trad_chinese":"渙","pinyin":"huàn","english":"Dispersing","binary":110010,"od":55},
    {"hex":60,"hex_font":"䷻","trad_chinese":"節","pinyin":"jié","english":"Restricting","binary":"010011","od":56},
    {"hex":61,"hex_font":"䷼","trad_chinese":"中孚","pinyin":"zhōngfú","english":"InnermostSincerity","binary":110011,"od":62},
    {"hex":62,"hex_font":"䷽","trad_chinese":"小過","pinyin":"xiǎoguò","english":"LittleExceeding","binary":"001100","od":61},
    {"hex":63,"hex_font":"䷾","trad_chinese":"既濟","pinyin":"jìjì","english":"AlreadyFulfilled","binary":"010101","od":64},
    {"hex":64,"hex_font":"䷿","trad_chinese":"未濟","pinyin":"wèijì","english":"NotYetFulfilled","binary":101010,"od":63}
]

from_binary_to_element_symbols = {
	'11': '🜁',
	'01': '🜂',
	'10': '🜄',
	'00': '🜃',
}

from_binary_to_element_ix = {
    '11': '-1',
    '01': '0',
    '10': '1',
    '00': '2',
}
