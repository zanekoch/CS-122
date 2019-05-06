import numpy as np

def edit_distance_matrix(s1, s2):
    #create matrix that is len(s1)xlen(s1)
    #fill with 0s
    outputMatrix = np.zeros((len(s1), len(s2)))

    for i in range(len(s1)):
        outputMatrix[i,0] = i
    for j in range(len(s2)):
        outputMatrix[0,j] = j

    for j in range(1, len(s2)): #i:0->len(s1)-1, s1 is vertical
        for i in range(1, len(s1)): #j:0->len(s1)-1, s2 is horizontal
            deletion = outputMatrix[i-1, j] + 1
            insertion = outputMatrix[i, j-1] + 1
            identity = outputMatrix[i-1, j-1] if s2[j] == s1[i] else np.inf
            substitution = outputMatrix[i-1, j-1] +1 if s2[j] != s1[i] else np.inf
            outputMatrix[i,j] = min(insertion,deletion, identity, substitution)
    return outputMatrix





if __name__ == '__main__':

    s1 = "HSAITYRKQKWHRMQIPKLCYYNLDERQTKNTTFGFRKTCHCCGCVGESCVTGCIIYWILRVFFMKQWSQNRKSFSLWGQLCSWQTCCKHSHPLAKAGDLIFSTQDPGAGCWQMTCDRWISTVSWQPQSIRRNPFTCQHVMFSKFNMIWARNMCVQTFICGLHGLYKLGLDMCVFRMGWRFIDCWNNSVTWYANSLTTARIFDTSSDASVWVMGSRHGRLEIPKKPRWFQQEDTFGGCPGEVFFFQFGAGAERAHRHLTCFWDWWFEAMQHWLSHPESCHKECIYKSPAQELFLNDIYKNIGHSEWHISPCVFGGCAVVRLCNSIEDDTFQFCSEGHGCNENAFPFIRDWEVTPAPGCHQIRACIAMRRATIKASHCKNIKRKIEQGHKELSQLCCQRVDGMSPKQKPNIHDYDPSASDYGFPAQEGMLGKEKMCYNQLDPYSMQNKQSYYGQIRNMVECFQEAWGMQTFTRCPWFEVCMTICMNFWTGNPADSYPLVTDILGCQSAYSVPRPCEPYRFWPEKSCIRKHIIWFFKPFQSAKPASDIWTDHPMGFHKCVWWSMQYTIIHQTRFELKWNKITMWTSLMQHKDNACHHCCTINEEEMGRGMNIPTQKIMGIFGDDLCGRFKIICKSAIDCMLWCDEQSPNYEPLACNPFFRHWKKHPCEQPIYTPQLNVKSTHAFRPITFETCHNCYMLTDQNCWEFGITDCVDDPKQWAMQAVCCHTSSCGGLYRINFYCHPETTDDDAQGMELSRHFNRVFSCSDRSYIGEHADKILQNQGFGVMLKNYLIHRDRQAGGYFHWTSIYFHGGSWNPPNGKDMMSHMLNFYQCLLSDQKSRGNVHPTYSLPASKMRTCNECMGYHFIPEEIWAMTPLVLLDQQCDYLIQVAHPFLRFPNSMYVVFDQFSHNDFSTTFKRETDDGLRLHSTRMRKNFMVPCLWWEPPGTLMNPSICSINYPAHLPLQITCTKPTEHLMCFYKHFMHVYIGKNRPLAKLEVDFQKSHTFYCYNFFYPMDWVNRWKGQCPGHAGYAKSGIFIQEKYQENWVADQFESYVNHGRKVPVLMVVPRKKEIYKCLCLWITQVAKTSYYWYRDFYKHPNFASWGCFHMSLKETNMDCERAHGRVDFLWEWGKSGCSEIREDRNDKKCNWDAEFSPHPCNTRHQVPLWPCLQFMWNDGAKAAEAPRHAEKIMMAFIKDWCSVCRINYPAYIFDCVVMHIEAREFDFIFKIKSQLGMVRWQRYVGFHTMNWKKEWSRCRTQCCNIDLEIGHPANRWREHAKDEYYRNYSQCDRAYTIQTFGRCQRACWQDSYFQDCVDPIWRRKSCCPMIMSREEINNHENKGRNMQACSWEMQTILVGNHFRSQRYHAEVIAGRVFFSCFWEIYAMTRGLTATLHFHFAEQDNRWIHPKRWLFIIMFYGQDVINQSVMYCNTTEYRWGLYQDENTVIRMGCDITCNQVDKDRQDMRQTWPFKTSSGITQSSTWAGQRCWEMDHHEYLSVQELYSGIIVLSMFHPDWNCNVEAQKNQMGTLLPSDFSCPPIADIVKIKSQKWMVPRACARRGWFGAPGNAWQRGEMMTWMLSRGIFRCNMSIISEWQYAFQCCPGLDMDLWVASQVDRKPGTETFWNHCHSGATKYNPRWPGKQMRRHVHIIHRRTSNMRGYKFNGQWYEADLMKCRAWSCQQKGAGLRASSIDHYTDYMLQNKPYSCHHCTTSLVSCMKFLFHEWMKTSCCSNCTNGPQCTPPDVMCDHKANEPFPWHVKLSRQIWPNTKPADDAMQLGKSHLQSIMWNMLPCCYRFTICFFMNAECNMTRGGGLFECPAIPPWREQKVALGRNCNCDKHYEKNQSAYKAVIPIVVIKSPFNFVFGMAIWWYMGAMNHQSNIQQMSYISSHRSDYWTFYDRTELGDFCWAMCWNARVRKHNRQEYKRVNSCCGKVTQAIMAMAWDPDSSRVVMFDGKFVWMAKFIGNGTKRQKCDWGYQDTTHVMYLCCACSSTFQKMFNPETLAHAYKCWHILFQLLKCQGMFLIEPHVSRARHIVEFDVIAMCLNVMLFQDSQYCEGVIVKMKDPWWCPKEKTKFAQCRTWNDQEHRYLWHAYGENEDWRFEKAPMRMVFESLDDSDQWQPHQKMGDASIDIFARGEQMSANDELFYWQPCAAEWFLKYKKHDRTIMETQYIHWNYCMGEIMDPVAMACHRHRPPKNALATAEMFTIPKYYKELEQYFLMRHMVDICWLDIIEASWTHCKPLKMESRKYHPRSWEISTGFITKLMLHIHRMNTLSEKAEAQHWCTIYRTHSPESFKSQWKDLSYWVNTNWQQAIDPSSVQERTGYPEMYNNFVMGIPFRHSPIANFTYTAKPQMAKEVFRTYGIKKWIWSECAHAFIGFREWFEMNISEQGNMPLWINCSLWNFAWNRWVCNITSYYALIHDGIEAPPDFQWSFYWQNTCLEFESGFHEQKIVKKALHYCAWTCKVESMNWPGYQEKVKHGVQHNMDFCGKPIKVWVHRQAYYRRDRYNEIPETEIYYCPWEMSRVIHDISIKETFGPFLDAPHRTFPDIMIRKKTDGTKLEPSIDGWPQQFQACLIFLNWDTSTFLLEYTENLTEEPVACPNYLKIAAACMGMHMGVNPMRCNRWWYLYRTPGHHGMKCDLVSFQFVGTMQCKNVETGGMFDPDPREDGMSWWHSIWNSFHNSPKYSSQPPDIFGTKDPRSVNMYVMPCMTLWHYKSLGKAVRQVCENTDDYEVVQRKRGQQVACSPVHRVVFKGAQSMTFTDKLPIEMSVWLGEFVDSCVCTSREDPISPCWIIQTTCGKKRMKLFQCWACSIPCDSFQYLRSEMHFNLRMNEACCQRIPVHWLFPEAQSYYWVGTQCWKIPIYMARTLKDEPCVFFECMGFYGAELLGFLMMSGNLFLGHVMGHQDWRVICDNGCRVAFQDVDPCCFNIIVWGSFAWKVSYTNEMMRYSDARVDADCWVMWSMSRQYSLADHVHQVGIMILADERSGCLWWPMIRRTQQEHCNNHIVICPAMVWPGILGQGWPLDQNLPMSCRVPNSYKPINSEQISYWNTPTQLTFWTCPSSDLFDYEGCPMCWSNGQEVDGQGLVYCTNDKHTPVTFFGVQLQQNCEFWVGAKMWRMTEYCALKFEEVEEVYPLDDIGLLRFLFDTRGLVQRRRRVKYVWANHYGECKHFMPCYFKEAWDYHARWENMSVIEQPQDVNEIGNKLRLHYRMVQGCWSHPLLPMDIRNTMKSVGIMSGPGNKRFNNYHHKETVGYNVEGHPCWEWHEAINNFHPRPVFSIWCCGYISARKFCDKKMNEPQDMQCIYLNWPYSLTYEAIDQAFMMDWGDYAWALNVMTDWKHISFRVLWMALCVFYANYWMYMWEICWQHIAKPAGCHIFQYARICDVQHSAVHCMHKSMYPVDEYETDSAYMPQIITLPEFNKTWHECGAHSRIAMAQVEMIEKGCRWTKQYTPFGIWMLHINWQCVYLSGQKDWLTHSDHRNQKVVQKSFYWPCYTPVKPIEQWEPYSEEPFQPYPLPNRPQQLMGQFFYFPAQHRTSGLIFTKDQEVAMCRIWNLHLHHWDQPTEMCWNTRWWDYWCCAYRAVFRWFYLLMFPPSDVMHPDLDCDMDPWVHLTEWDWIFEILWILELCCSNTMFNQRMVIEGSFMHSPSELYDRGQDPQYPSSSRYATCFYPYAYLKYAHSKFQKYWLIGGLKTFSIMTNYVMTIHHSGPIACLGLQYRNWCQKAQESTFKSYRPTGNWCGNDKRGTNHMDIFNVRGFAHMGWTVMTVISEGTECGPCVQCTLCELCERYAQDDYMLRYKHECQPMFGTCHMWCELDCTHMDNSGVSAMMIKISEQGSVYSWICAVKCFARIVSIRASGGKDNLYPYNGPKNSGTSEMFGSLLIIQSVMAHQQLAPQEPVENTRGRQHVNPSEMEANNNTSTICYYYMSKHVPLETFYAHEVFESGRYFRRADKGRLNLNGVVSEFEQEFQSTFQCHYVRMDRDGLWWQDAQVCGPSARCPWGDYHDSVSYEWGFHDPTIWNGGWFKTHTKKRGSACSIDKNLLMHKTQDEWKFYCTNGLQHCDNDWYYFYVKPSLTHYIQTSQSDVTKIIRTPNNCSLYYYGSYQHGLDTDYLKRANGHVIVVWDICDFHPFECLEWALAAACNTVVICVWPDIRNIWVCMKNWLGRCWPMGDYDYATGSGNMEQYSCSWTEQKRFQPGCQQPHQKPASSQLPDSWVCEHKPATCVYAIMCPIGSAKKCMCQWYQQNQHCAWAIPDAEAIMEGWLVWTYGRLYMDFELIMYPCLVMDANACMVASMCINWNRKTDAWIAIQYFPRKRDTAHEDYVEGCGKCGAININTRAVNDMKEKKHTWPEENMTRMQEKVRPA"
    s2 = "HSAITYCPLKQKWHRMQIPKLCYYNLIAERQTKNTTFGAGRKTCICCGCVRIAYGHTWHSCVTGCIIYLILRVFFMKQWSQNIKSASLWSQLCSWQPCCKHSHPLLCSAFVWMAIGDLIDSTQDPGAGCATQHTCDRWIKTVVSHPQDQRFKEDYIMNMTCQHVMFSKANMIWARNMTFDFSLYKLGLDMCVRDGWRIDCWNNSAWLTQWYAREFDTSMGSRAQCLMQMTGRNNMKEIPKKPRWFQCEDTFGGCPGEVFFFQFGAGAQCRMWLKIQFYWLSHPESCHKEQIYKSPAQKPKYKNIGHSQVQTRWSPCVFGGCAKFLNLLCNSISGFCCNENAFPFIRDAEVTPALGCHQIRACIAMRRATIANDWWAECASHCKYIKRKIRDQGHAGAHENLSQLCCCRVDGMFPKQKPIIHDHDPSAPAQEGILGKQKMCFFQLDPYNMQNIQSYYGQRNCKCPDHPMCIGRICMNIPCIQTGNAADSWCSVQTPPLVRDIEGCQSAYSIPRPCEPYRFWPEKRKHQMCHDFQIKIWFFKPFQSAKPVCEFENGPQRRCWPCDMGLHKCVTIIHQTRFTLKWNKITMWTSKAHLFCPLDQPGHCSIRWRYTTMNEEEMGRGIMGIFADDLCMRFLCRIICKSEIYCMLWCDEQSPNYEPLACNPFFRMWKLHDVKSPITFENIMYPHNCYMLTWQNDWEFGDTDCVGDPKQWLMQARLCCHTSSCGGLYRINFYCHGMELSRYFHRVFSCSDRSYIGEGLMKILQNQFMFVVLKIFIEPCALGYLIHRDRQAGGYFCITVRYKWTSIYFHPNGKDMMSHLRDSQKSRLNVWTWKWMCPVTYSLPSKHRRFCNECMGYHFIPEEIWAMTPKVLLDQQCLYLHPFLRFPNGMYDPQIQVRDQLSINDFSTGFKRETDDGLRLHSTRMRKNHNAPFHMVGTLDMCFYIASGPHLHGPHKICSIRYPAHLPLQITKTKMTEHLIYHSTNFMHVFSIIHPEDDIGKNRPLAKLEVDFQKSHEFRPLEWYPMDWVWQKRTARSHAGYAKSGIFIQEKYQENWDADGHCKDANSYVNHGRHVPHIVPWKRCGMVVPRKKCLCLWITQVAKTSYYWYRDFYKHVPKGVAEFFASWGCIQVVDHMSIKETNRATGNRAGTFDFLWEWGKFWEGCFVPGCSQVNAWPCNNDAEFSPHPCNTRHQVPLWPCLQFMWAAEAPAEKIMMAFIKDRINYPAMIKSQVRSKDCVVMHIEAREFDFIFRWRIIPIKSQLGSVRWQRAVGFHGENWKKSQCWTKHSRSRCRTQGWWVCWGTQAIGIDLEWYIAQFFREGHITNRWREHAKDEYYQILAITEQTFGRCQRAEYFQNMYKLGGRKSTQDSGTGCPMMMSREDWNNNMQAQLHSWEMQTILVREHAVFRSQRYHAEVISFNRSFFSCFWEIYAMTRGLTANRWMHFFIAEMCDVVNQSTEKLTFQDKRHRKYYRWGLTVIPHQQWNFIRRWTSMGCDEWPRDCHKTCNQVDRQLVLYPKHNMRQTWPFKTSSGITQSSTWAGHRCWEMDHHEYEVQELYTEIHKMIMVPWDNCNVEAQKNQMGTLLPSDFSCPPIVRIKSQKWMVPRACARRHWFGAPGNAWQRRMMTWMLSRNGSCVRQTSYAFIIHQCCPGLDMDLWLASQVDRKPGTETHSGATKYKPRWPKKQMHVYIIIRRTSNMRMINEDWACYKHNGQWYKCVAWSIFIQTVTQQKGAGLRASSIDHYMLQNSCHSCTTSIVRCMKFPFHEMKMSCCSNCTNGMQGQDVMWVGYHQRLDWANEPFPWHPKLSVQIWPNTKPADDAMQLGKLHLQSIMPCCYRFTICFFMNAELNMTRGGQCLFECPTIHPWREHKVAKGRNCQTSSWSKDCDKHYECYFTCNDTEGHCAKAHKAHYPIVVINFVFNMAYMGAMNHRTESNISSHRSDHWTKDGTEMGDFCWAIQVQVCWNARVRKHNRDGRANAVPRSRTNLYCCGFAGWQCYVTQAIMAMASHCYDPCSSFRSNMRDFEVPMFDGKFVTKRQKCDWGYQDVKCLPDNGYECKDTGDQACSSTMKIYAIKQHMFNIEPLAHAYKCWYMELFQCQGMWLIEPHVSRKIARHIMCLNVMLEQDSQYCEGVIVKMKDPEKTKQAQCRVWPRPNTYEWNDLMHRYLWHHYVNFSPFYGENEDWRFEKAPMRMVFESLDSDQWPGDQSARGEKLRMSANDCIMHHGMVVLFYWQICFAIAHWFLKYKKHDITHMETQYIHWNYCMGEIMDPVETIPQYFMMRHMVDICWLDIIEASWTHCKPLKMESRKYHPRSWEISTNNFITKLMLHIHRMMCDHCYPQTLAQHWCAYFAKQHNIYRTQMKHTKSPESFKSQWKDPVNTNWQQRSIDPSSVQERTGYPEMYNNQVMGIPRISPNPYEKQVEFFRTYGIKKWIWSAHAGFVIRDANISEQGHMTLWTNCSLWNFAWNWWQWAMNYRATSYYVSDKLQDQLIFDGIPDFQTPPDFQWSFYWQNTCLEFESGFHEQWRWWQRIVKKALPSMNMPGYQAKVKHGVQHGKNHKVWVHRCAYYRRDRYTMFWFYIAIPETEIYYCPWEMMLNGCRGISDYSIKETFGPFLDAPFTFPDIMIRKHPSRMGTKLEPSIDGWQFQACLIFLNTKTELLMYTPQGNLAEDPVACHYNYPKIAAACMGMCMGVNPMRQNRWWLLYSFQMQCKNVAECQHGGMFDPIPREDGMSKWHSIRFNSFHNSPKGSSTNQSQWFGWKNPRQKWDGSMVNKIQHYVMPWHYKILGKQVCENTEEQRNKNECDNVVQNKRGQQVACETVRFVHRVVEKGAQSMTFTDKYPAAMSVWLWRFVDSKCFYWFYVCTSRETEQISYPISPRANRIEEMMKIIQTTCGKKRMKLFQCLLTSQCYKACSDHMTMPCDSFTYLEDVTYLRSEMHFNLRMNEPMKACWAHWFFPEANGSYYNVGTQCWHIPIYMARTLKDEFECMGFYGAELLGFCMMSGNTKVWHQYQFLDWRVICDNGCRVAFQDVDPCCFNISVWGSIESPECTAWKVSYTNEMPAYSDARVDADCWVMWSHSRQYSLADHLADERSGCLWWPLIRRTQQEHCNNHIVICPAMVWPGILFYSVANLPSCRVPNSYKPIHNLMWWYIVQITFWTCLSSDVIGFIVSQTIDEYHIPKMLLLNIYNCWSNGFEVDGQGGGGVYCTNDKHTPVTFFGVQLQQNCEFIVGAMMWRETPLDMLDMSPYCVAIELRFLRFDFGQRRCRVKYVWALHCGECKHFMPCYSMKEAWDIHARWENPNISGYTCFKAKQGMQPQDVGKLRLHARMVQGCWSHPTLKGNKSVGIMSGPGNKRFNNYAYENEFVHKETVGYNEWHEAINNRHPRPNNNMKKGKWISARPIWKESCDKKMNEPQDMQCIYLNWEAIDQAFMWGDYAWPLNVMTDRVLWMFYANMLPANRWMYMWEICWQHRAKPANCHIFFNRYARICDVQEFVCGKKHSFDIMMIRWEKDHEDECMHKSMYPVDEYDDTEGDSAYMEFAHSRIAMAQVEMIEKGCRMINWQCVTIYTFLSGQKDWLTSSKLSFMHNQKVLPMCKQKHYFYHPCYTEQWEPYSEERPFPAQHRTFIFTLTHWWWMWNTKDQEVHMCRIWNLHLHYWDQPTEMCRHLDYECWCGAYRAHFSWFYLLMDCACDFDPWVHLTEWDELIEMFNQRMVIEEKFFDMFMHSPSEDRGQDPQYPSSSFPYCIKMDQMAYLKYATSTLQKYWLIGGLKTFCMCMYWIMTNYVMTIHHSKWRPMHDNEFEQNPEVYRHYNWCQKGQESTKDQCHKWECCKAQGNGKISFNVRGFAHMGWTVMTVEGTECGPCCNMYAGQLTLCELCERCKQDDYMQGDQYDPECQPMFGTCWQTHMDNSGNSQKMMAIDIMMIKVYSWKCIRSFGKPKNTSEMFGSLLIILAPEPVENTSGRQHCYYYMSKHVPLRCLGTKYAIEIQKWKVNEYFRRAFKGRLNLNGVVSEFEQEFISLQFQCHYVRMERDGLCWVCGPSARCPWSIYHDTITSYEWGFHMQTIWNGGKKNLLMKYSVKWKTNGAKNDFYVKPSLFHSIQWSSDVFKIIRTPNNCSLTVVIYGSYTHGLDTDYLKRANGHVIVVMPLEDCHECLEKALAAACNTVVICVWPDIRNIWMCMKNCLWDYDYSQMVPPHDWMSGNMEQYSCSWTEQKPGCQRCPHQDPAFSQLPDIVWVCEHCPNTCVYAIMCPFGEVNTLMHQWYQQAHCAWAICDAEAIMEGWLVWTYGRLYMRFELIVMDANASGHSHKKTDAWCHIYKAIQYEWRKRDTAHEDGVEGCGKCGASRTKMWFINTRASNPVNNMQNDLQMGKEKKHTWPEENMTRHQEKVMSRVA"
    output = edit_distance_matrix(s1, s2)
    print(output)
    print(output[len(s1)-1, len(s2)-1])