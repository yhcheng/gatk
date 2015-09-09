package org.broadinstitute.hellbender.utils;

import htsjdk.variant.vcf.VCFHeader;

import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.utils.SampleUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collection;
import java.util.HashMap;
import java.util.Set;


/**
 * Testing framework for sample utilities class.
 *
 */
public class SampleUtilsUnitTest extends BaseTest {
    @Test(expectedExceptions=UserException.class)
    public void testBadSampleFiles() throws Exception {
        HashMap<String, VCFHeader> headers = new HashMap();
        //GATKVariantContextUtils.GenotypeMergeType);
        // create headers
        Set<String> samples = SampleUtils.getSampleList(headers, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);
    }
}
