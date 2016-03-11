package org.broadinstitute.hellbender.utils;

import htsjdk.tribble.index.Index;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class IndexUtilsUnitTest extends BaseTest{

    @DataProvider(name= "okFeatureFiles")
    public Object[][] okFeatureFiles() {
        return new Object[][] {
                { new File(getToolTestDataDir(), "test_variants_for_index.vcf")},
                { new File(getToolTestDataDir(), "test_variants_for_index.gvcf.vcf")},
                { new File(getToolTestDataDir(), "test_bed_for_index.bed")},
        };
    }

    @Test(dataProvider = "okFeatureFiles")
    public void testLoadIndex(final File featureFile) throws Exception {
        final Index index = IndexUtils.loadTribbleIndex(featureFile);
        Assert.assertNotNull(index);
    }

    @Test(expectedExceptions = UserException.class)
    public void testLoadIndexFail_tooNew() throws Exception {
        final File featureFile = new File(getToolTestDataDir(), "test_variants_for_index.newerThanIndex.vcf");
        featureFile.setLastModified(System.currentTimeMillis()); //touch the file but the index
        final Index index = IndexUtils.loadTribbleIndex(featureFile);
        Assert.assertNotNull(index);
    }

    @Test
    public void testLoadIndexFail_noIndex() throws Exception {
        final File featureFile = new File(getToolTestDataDir(), "test_variants_for_index.noIndex.vcf");
        final Index index = IndexUtils.loadTribbleIndex(featureFile);
        Assert.assertNull(index);
    }

}
