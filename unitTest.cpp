#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <string>
#include <boost/test/unit_test.hpp>

#include "vcfData.hpp"
#include "vcfParser.hpp"

using namespace vcf;

BOOST_AUTO_TEST_CASE( testSplitWithDelimiter )
{
    VcfParser parser("");
    auto res1 = parser.splitWithDelimiter("hello, world", ',');
    BOOST_CHECK_EQUAL(res1.size(), 2);
    BOOST_CHECK_EQUAL(res1[0].compare("hello"), 0);
    BOOST_CHECK(res1[0].compare(" hello") != 0);
    BOOST_CHECK_EQUAL(res1[1].compare("world"), 0);

    auto res2 = parser.splitWithDelimiter("hello", ',');
    BOOST_CHECK_EQUAL(res2.size(), 1);
    BOOST_CHECK_EQUAL(res2[0].compare("hello"), 0);
    
    auto res3 = parser.splitWithDelimiter("", ',');
    BOOST_CHECK(res3.empty());
    
    auto res4 = parser.splitWithDelimiter("hello, world, how, are you?", ',');
    BOOST_CHECK_EQUAL(res4.size(), 4);
}


BOOST_AUTO_TEST_CASE( testSplitAndCapitalise )
{
    VcfParser parser("");
    auto res1 = parser.splitAndCapitalise("ATt, GcT");
    BOOST_CHECK_EQUAL(res1.size(), 2);
    BOOST_CHECK_EQUAL(res1[0].compare("ATT"), 0);
    BOOST_CHECK(res1[0].compare("ATt") != 0);
    BOOST_CHECK_EQUAL(res1[1].compare("GCT"), 0);

    auto res2 = parser.splitAndCapitalise("ATt3df, jklo, p;Y");
    BOOST_CHECK_EQUAL(res2.size(), 3);
    BOOST_CHECK_EQUAL(res2[0].compare("ATT3DF"), 0);
    BOOST_CHECK_EQUAL(res2[1].compare("JKLO"), 0);
    BOOST_CHECK_EQUAL(res2[2].compare("P;Y"), 0);
    
    auto res3 = parser.splitAndCapitalise("");
    BOOST_CHECK(res3.empty());
}

BOOST_AUTO_TEST_CASE( testExtractGeneName )
{
    VcfParser parser("");
    std::string str1 = "TCR=22;TR=44;WE=68757;WS=68739;Gene=ENSG00000178591;SYMBOL=DEFB125;CSN=1;BIOTYPE=protein_coding";
    auto res1 = parser.extractGeneName(str1);
    BOOST_CHECK_EQUAL( res1.compare("ENSG00000178591"), 0 );
    BOOST_CHECK( res1.compare("") != 0 );
    
    auto res2 = parser.extractGeneName("WS=68739;aGene=ENSG00000178591;SYMBOL=DEFB125;");
    BOOST_CHECK_EQUAL( res2.compare(""), 0 );
    
    auto res3 = parser.extractGeneName("WS=68739;Gene=;SYMBOL=DEFB125;");
    BOOST_CHECK_EQUAL( res3.compare(""), 0 );
    
    auto res4 = parser.extractGeneName("Gene=ENSG00000178591;SYMBOL=DEFB125;");
    BOOST_CHECK_EQUAL( res4.compare("ENSG00000178591"), 0 );
    
}

BOOST_AUTO_TEST_CASE( testVcfParser )
{
    VcfParser parser("testData/testVcfFile.vcf");
    BOOST_CHECK(parser.hasNextValidRecord());
    
    auto validRecord1 = parser.getNextValidRecord();
    BOOST_CHECK_EQUAL( validRecord1.chrom.compare("20"), 0 );
    BOOST_CHECK_EQUAL( validRecord1.pos, 65900 );
    BOOST_CHECK_EQUAL( validRecord1.id.compare("rs6053810"), 0 );
    BOOST_CHECK_EQUAL( validRecord1.ref.size(), 1 );
    BOOST_CHECK_EQUAL( validRecord1.ref[0].compare("G"), 0 );
    BOOST_CHECK_EQUAL( validRecord1.alt.size(), 1 );
    BOOST_CHECK_EQUAL( validRecord1.alt[0].compare("A"), 0 );
    BOOST_CHECK( validRecord1.pass );
    BOOST_CHECK_EQUAL( validRecord1.geneName.compare("ENSG00000178591"), 0 );
    
    BOOST_CHECK(parser.hasNextValidRecord());
    auto validRecord2 = parser.getNextValidRecord();
    BOOST_CHECK_EQUAL( validRecord2.pos, 67500 );   // has skipped the invalid record
    BOOST_CHECK_EQUAL( validRecord2.alt.size(), 1 );    // and parsed correctly the long alt field
    BOOST_CHECK_EQUAL( validRecord2.alt[0].compare("TTGGTATCTAG"), 0 );
    
    BOOST_CHECK(parser.hasNextValidRecord());
    auto validRecord3 = parser.getNextValidRecord();
    BOOST_CHECK_EQUAL( validRecord3.alt.size(), 2 );
    BOOST_CHECK_EQUAL( validRecord3.alt[0].compare("C"), 0 );
    BOOST_CHECK_EQUAL( validRecord3.alt[1].compare("A"), 0 );

    BOOST_CHECK(parser.hasNextValidRecord());
    BOOST_CHECK(parser.hasNextValidRecord()); // can call multiple times, will not change until data requested
    auto validRecord4 = parser.getNextValidRecord();
    
    BOOST_CHECK( !parser.hasNextValidRecord() );
}

