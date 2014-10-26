#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <string>
#include <boost/test/unit_test.hpp>

#include "vcfData.hpp"
#include "vcfParser.hpp"

using namespace vcf;

BOOST_AUTO_TEST_CASE( testSplitWithDelimiter )
{
    VcfParser parser("testData/testVcfFile.vcf");
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

BOOST_AUTO_TEST_CASE( testExtractGeneName )
{
    VcfParser parser("testData/testVcfFile.vcf");
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
    BOOST_CHECK_EQUAL( validRecord1.ref.compare("G"), 0 );
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

BOOST_AUTO_TEST_CASE( testIsIncluded )
{
    BOOST_CHECK(GeneData::isIncluded("", ""));
    BOOST_CHECK(GeneData::isIncluded("a", "a"));
    BOOST_CHECK(!GeneData::isIncluded("a", "A"));
    BOOST_CHECK(!GeneData::isIncluded("a", "ab"));
    BOOST_CHECK(GeneData::isIncluded("ab", "a"));
    BOOST_CHECK(GeneData::isIncluded("ABCD", "AB"));
    BOOST_CHECK(GeneData::isIncluded("ABCD", "BC"));
    BOOST_CHECK(GeneData::isIncluded("ABCD", "CD"));
    BOOST_CHECK(!GeneData::isIncluded("ABCD", "CDE"));
    BOOST_CHECK(!GeneData::isIncluded("ABCD", "BD"));
    BOOST_CHECK(GeneData::isIncluded("ABCD", ""));
}

BOOST_AUTO_TEST_CASE( testEvaluateMutationType)
{
    BOOST_CHECK(GeneData::evaluateMutationType("", "") == MutationType::IDEN);
    BOOST_CHECK(GeneData::evaluateMutationType("A", "A") == MutationType::IDEN);
    BOOST_CHECK(GeneData::evaluateMutationType("AT", "AT") == MutationType::IDEN);
    BOOST_CHECK(GeneData::evaluateMutationType("A", "G") == MutationType::SVN);
    BOOST_CHECK(GeneData::evaluateMutationType("C", "T") == MutationType::SVN);
    BOOST_CHECK(GeneData::evaluateMutationType("AT", "A") == MutationType::DEL);
    BOOST_CHECK(GeneData::evaluateMutationType("ABGHG", "H") == MutationType::DEL);
    BOOST_CHECK(GeneData::evaluateMutationType("ABHG", "ABH") == MutationType::DEL);
    BOOST_CHECK(GeneData::evaluateMutationType("BHSJ", "HSJ") == MutationType::DEL);
    BOOST_CHECK(GeneData::evaluateMutationType("D", "BDHA") == MutationType::INS);
    BOOST_CHECK(GeneData::evaluateMutationType("D", "DNJ") == MutationType::INS);
    BOOST_CHECK(GeneData::evaluateMutationType("DB", "ADBK") == MutationType::INS);
    BOOST_CHECK(GeneData::evaluateMutationType("DB", "DBHJ") == MutationType::INS);
    BOOST_CHECK(GeneData::evaluateMutationType("DB", "PLDB") == MutationType::INS);
    BOOST_CHECK(GeneData::evaluateMutationType("DB", "") == MutationType::DEL);
    BOOST_CHECK(GeneData::evaluateMutationType("DB", "ADKBHJ") == MutationType::MVN);
    BOOST_CHECK(GeneData::evaluateMutationType("DB", "ADKBHJ") == MutationType::MVN);
    BOOST_CHECK(GeneData::evaluateMutationType("AAAA", "ATAT") == MutationType::MVN);
    BOOST_CHECK(GeneData::evaluateMutationType("CTCT", "CACTTTTTT") == MutationType::MVN);
}

BOOST_AUTO_TEST_CASE( testGenomeData )
{
    GenomeData gd("testData/testVcfFile2.vcf");

    // all records get correctly parsed and aportioned to genes
    BOOST_CHECK_EQUAL(gd.m_positionRecords.size(), 6);  // 8 in total, 2 do not PASS
    
    // check genes
    BOOST_CHECK_EQUAL(gd.m_geneData.size(), 2);     // 2 genes
    BOOST_CHECK(gd.m_geneData.count("ENSG00000178591"));
    BOOST_CHECK(gd.m_geneData.count("ENSG00000125788"));
    auto geneMutCounter1 = gd.m_geneData.at("ENSG00000178591").geneMutationTypeCounter().counter;
    auto geneMutCounter2 = gd.m_geneData.at("ENSG00000125788").geneMutationTypeCounter().counter;
    BOOST_CHECK_EQUAL(geneMutCounter1[MutationType::SVN], 3);
    BOOST_CHECK_EQUAL(geneMutCounter1[MutationType::INS], 1);
    BOOST_CHECK_EQUAL(geneMutCounter1[MutationType::DEL], 1);
    BOOST_CHECK_EQUAL(geneMutCounter1[MutationType::MVN], 0);
    BOOST_CHECK_EQUAL(geneMutCounter1[MutationType::IDEN], 0);

    BOOST_CHECK_EQUAL(geneMutCounter2[MutationType::SVN], 0);
    BOOST_CHECK_EQUAL(geneMutCounter2[MutationType::INS], 0);
    BOOST_CHECK_EQUAL(geneMutCounter2[MutationType::DEL], 2);
    BOOST_CHECK_EQUAL(geneMutCounter2[MutationType::MVN], 1);
    BOOST_CHECK_EQUAL(geneMutCounter2[MutationType::IDEN], 0);
    
    // check totals
    auto totalCounter = gd.m_mutationCounter.counter;
    BOOST_CHECK_EQUAL(totalCounter[MutationType::SVN], 3);
    BOOST_CHECK_EQUAL(totalCounter[MutationType::INS], 1);
    BOOST_CHECK_EQUAL(totalCounter[MutationType::DEL], 3);
    BOOST_CHECK_EQUAL(totalCounter[MutationType::MVN], 1);
    BOOST_CHECK_EQUAL(totalCounter[MutationType::IDEN], 0);
}
