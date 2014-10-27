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

BOOST_AUTO_TEST_CASE( testExtractGeneNames )
{
    VcfParser parser("testData/testVcfFile.vcf");
    std::string str1 = "TCR=22;TR=44;WE=68757;WS=68739;Gene=ENSG00000178591;SYMBOL=DEFB125;CSN=1;BIOTYPE=protein_coding";
    auto res1 = parser.extractGeneNames(str1);
    BOOST_CHECK_EQUAL( res1.size(), 1 );
    BOOST_CHECK_EQUAL( res1.at(0).compare("ENSG00000178591"), 0 );
    BOOST_CHECK( res1.at(0).compare("") != 0 );
    
    auto res2 = parser.extractGeneNames("WS=68739;aGene=ENSG00000178591;SYMBOL=DEFB125;");
    BOOST_CHECK( res2.empty() );
    
    auto res3 = parser.extractGeneNames("WS=68739;Gene=;SYMBOL=DEFB125;");
    BOOST_CHECK( res3.empty());
    
    auto res4 = parser.extractGeneNames("Gene=ENSG00000178591,ENSG00000222,ENSG000001111;SYMBOL=DEFB125;");
    BOOST_CHECK_EQUAL( res4.size(), 3 );
    BOOST_CHECK_EQUAL( res4.at(0).compare("ENSG00000178591"), 0 );
    BOOST_CHECK_EQUAL( res4.at(1).compare("ENSG00000222"), 0 );
    BOOST_CHECK_EQUAL( res4.at(2).compare("ENSG000001111"), 0 );
    
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
    BOOST_CHECK_EQUAL( validRecord1.geneNames.size(), 1 );
    BOOST_CHECK_EQUAL( validRecord1.geneNames.at(0).compare("ENSG00000178591"), 0 );
    
    BOOST_CHECK(parser.hasNextValidRecord());
    auto validRecord2 = parser.getNextValidRecord();
    BOOST_CHECK_EQUAL( validRecord2.pos, 67500 );   // has skipped the invalid record
    BOOST_CHECK_EQUAL( validRecord2.alt.size(), 1 );    // and parsed correctly the long alt field
    BOOST_CHECK_EQUAL( validRecord2.alt[0].compare("TTGGTATCTAG"), 0 );
    BOOST_CHECK_EQUAL( validRecord2.geneNames.size(), 3 );
    BOOST_CHECK_EQUAL( validRecord2.geneNames.at(0).compare("ENSG00000178591"), 0 );
    BOOST_CHECK_EQUAL( validRecord2.geneNames.at(1).compare("PurpleEyes"), 0 );
    BOOST_CHECK_EQUAL( validRecord2.geneNames.at(2).compare("GreenHair"), 0 );

    
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
    BOOST_CHECK(PositionRecordProcessor::isIncluded("", ""));
    BOOST_CHECK(PositionRecordProcessor::isIncluded("a", "a"));
    BOOST_CHECK(!PositionRecordProcessor::isIncluded("a", "A"));
    BOOST_CHECK(!PositionRecordProcessor::isIncluded("a", "ab"));
    BOOST_CHECK(PositionRecordProcessor::isIncluded("ab", "a"));
    BOOST_CHECK(PositionRecordProcessor::isIncluded("ABCD", "AB"));
    BOOST_CHECK(PositionRecordProcessor::isIncluded("ABCD", "BC"));
    BOOST_CHECK(PositionRecordProcessor::isIncluded("ABCD", "CD"));
    BOOST_CHECK(!PositionRecordProcessor::isIncluded("ABCD", "CDE"));
    BOOST_CHECK(!PositionRecordProcessor::isIncluded("ABCD", "BD"));
    BOOST_CHECK(PositionRecordProcessor::isIncluded("ABCD", ""));
}

BOOST_AUTO_TEST_CASE( testEvaluateMutationType)
{
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("", "") == MutationType::IDEN);
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("A", "A") == MutationType::IDEN);
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("AT", "AT") == MutationType::IDEN);
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("A", "G") == MutationType::SVN);
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("C", "T") == MutationType::SVN);
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("AT", "A") == MutationType::DEL);
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("ABGHG", "H") == MutationType::DEL);
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("ABHG", "ABH") == MutationType::DEL);
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("BHSJ", "HSJ") == MutationType::DEL);
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("D", "BDHA") == MutationType::INS);
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("D", "DNJ") == MutationType::INS);
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("DB", "ADBK") == MutationType::INS);
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("DB", "DBHJ") == MutationType::INS);
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("DB", "PLDB") == MutationType::INS);
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("DB", "") == MutationType::DEL);
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("DB", "ADKBHJ") == MutationType::MVN);
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("DB", "ADKBHJ") == MutationType::MVN);
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("AAAA", "ATAT") == MutationType::MVN);
    BOOST_CHECK(PositionRecordProcessor::evaluateMutationType("CTCT", "CACTTTTTT") == MutationType::MVN);
}

BOOST_AUTO_TEST_CASE( testGenomeData )
{
    GenomeData gd("testData/testVcfFile2.vcf");
    
    BOOST_CHECK_EQUAL(gd.m_geneMutationCounter->size(),3);
    BOOST_CHECK(gd.m_geneMutationCounter->count("ENSG00000178591"));
    BOOST_CHECK(gd.m_geneMutationCounter->count("ENSG00000125788"));
    BOOST_CHECK(gd.m_geneMutationCounter->count("BlueEyes"));
    auto counter1 = gd.m_geneMutationCounter->at("ENSG00000178591").counter;
    auto counter2 = gd.m_geneMutationCounter->at("ENSG00000125788").counter;
    auto counter3 = gd.m_geneMutationCounter->at("BlueEyes").counter;
    BOOST_CHECK_EQUAL(counter1[MutationType::SVN], 5);
    BOOST_CHECK_EQUAL(counter1[MutationType::INS], 2);
    BOOST_CHECK_EQUAL(counter1[MutationType::DEL], 2);
    BOOST_CHECK_EQUAL(counter1[MutationType::MVN], 1);
    BOOST_CHECK_EQUAL(counter1[MutationType::IDEN], 0);
    
    BOOST_CHECK_EQUAL(counter2[MutationType::SVN], 2);
    BOOST_CHECK_EQUAL(counter2[MutationType::INS], 1);
    BOOST_CHECK_EQUAL(counter2[MutationType::DEL], 2);
    BOOST_CHECK_EQUAL(counter2[MutationType::MVN], 2);
    BOOST_CHECK_EQUAL(counter2[MutationType::IDEN], 0);
    
    BOOST_CHECK_EQUAL(counter3[MutationType::SVN], 2);
    BOOST_CHECK_EQUAL(counter3[MutationType::INS], 1);
    BOOST_CHECK_EQUAL(counter3[MutationType::DEL], 1);
    BOOST_CHECK_EQUAL(counter3[MutationType::MVN], 0);
    BOOST_CHECK_EQUAL(counter3[MutationType::IDEN], 0);
    
    BOOST_CHECK_EQUAL(gd.m_totalMutationCounter->totalNumberOfMutations(), 14);
    auto totalCounter = gd.m_totalMutationCounter->counter;
    BOOST_CHECK_EQUAL(totalCounter[MutationType::SVN], 6);
    BOOST_CHECK_EQUAL(totalCounter[MutationType::INS], 2);
    BOOST_CHECK_EQUAL(totalCounter[MutationType::DEL], 4);
    BOOST_CHECK_EQUAL(totalCounter[MutationType::MVN], 2);
    BOOST_CHECK_EQUAL(totalCounter[MutationType::IDEN], 0);
}

