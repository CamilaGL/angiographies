'''Evaluation script to get metrics óf segmentations against ground truth labels'''


import numpy as np
import os
import collections
import inspect
import json
import hashlib
from datetime import datetime
from multiprocessing.pool import Pool
import numpy as np
#import pandas as pd
import SimpleITK as sitk
from angiographies.utils.iojson import writejson
from collections import OrderedDict

# ----------------------------------------------- metrics

def assert_shape(test, reference):

    assert test.shape == reference.shape, "Shape mismatch: {} and {}".format(
        test.shape, reference.shape)


class ConfusionMatrix:

    def __init__(self, test=None, reference=None):

        self.tp = None
        self.fp = None
        self.tn = None
        self.fn = None
        self.size = None
        self.reference_empty = None
        self.reference_full = None
        self.test_empty = None
        self.test_full = None
        self.set_reference(reference)
        self.set_test(test)
        #self.skirefcl = None
        #self.skitestcl = None

    def set_test(self, test):

        self.test = test
        self.reset()

    def set_reference(self, reference):

        self.reference = reference
        self.reset()

    # ------ added by Camila
    def set_testcl(self, testcl):

        self.testcl = testcl
        self.resetcl()

    def set_referencecl(self, referencecl):

        self.referencecl = referencecl
        self.resetcl()

    def resetcl(self):
        self.clp2vollintersect = None
        self.clp2volltotalcl = None
        self.cll2volpintersect = None
        self.cll2volptotalcl = None

    def set_skirefcl(self, refcl):
        self.skirefcl=refcl
        self.resetclski()
    
    def set_skitestcl(self, testcl):
        self.skitestcl=testcl
        self.resetclski()

    def resetclski(self):
        self.clp2vollintersectski = None
        self.clp2volltotalclski = None
        self.cll2volpintersectski = None
        self.cll2volptotalclski = None

    # def computeSkiCl(self):
    #     if self.test is not None:
    #         self.skitestcl = skeletonize(self.test)
    #         if self.skitestcl is None:
    #             print("test no skeletonizo")
    #     if self.reference is not None:
    #         self.skirefcl = skeletonize(self.reference)
    #         if self.skirefcl is None:
    #             print("ref no skeletonizo")
    #     if (self.skirefcl is None):
    #         print("ref no está")
    #     if (self.skitestcl is None):
    #         print("test no está")



    # ------------------------
    

    def reset(self):

        self.tp = None
        self.fp = None
        self.tn = None
        self.fn = None
        self.size = None
        self.test_empty = None
        self.test_full = None
        self.reference_empty = None
        self.reference_full = None

    def compute(self):

        if self.test is None or self.reference is None:
            raise ValueError("'test' and 'reference' must both be set to compute confusion matrix.")

        assert_shape(self.test, self.reference)

        self.tp = int(((self.test != 0) * (self.reference != 0)).sum())
        self.fp = int(((self.test != 0) * (self.reference == 0)).sum())
        self.tn = int(((self.test == 0) * (self.reference == 0)).sum())
        self.fn = int(((self.test == 0) * (self.reference != 0)).sum())
        self.size = int(np.prod(self.reference.shape, dtype=np.int64))
        self.test_empty = not np.any(self.test)
        self.test_full = np.all(self.test)
        self.reference_empty = not np.any(self.reference)
        self.reference_full = np.all(self.reference)

    def get_matrix(self):

        for entry in (self.tp, self.fp, self.tn, self.fn):
            if entry is None:
                self.compute()
                break

        return self.tp, self.fp, self.tn, self.fn

    def get_size(self):

        if self.size is None:
            self.compute()
        return self.size

    def get_existence(self):

        for case in (self.test_empty, self.test_full, self.reference_empty, self.reference_full):
            if case is None:
                self.compute()
                break

        return self.test_empty, self.test_full, self.reference_empty, self.reference_full


def dice(test=None, reference=None, confusion_matrix=None, nan_for_nonexisting=True, **kwargs):
    """2TP / (2TP + FP + FN)"""

    if confusion_matrix is None:
        confusion_matrix = ConfusionMatrix(test, reference)

    tp, fp, tn, fn = confusion_matrix.get_matrix()
    test_empty, test_full, reference_empty, reference_full = confusion_matrix.get_existence()

    if test_empty and reference_empty:
        if nan_for_nonexisting:
            return float("NaN")
        else:
            return 0.

    return float(2. * tp / (2 * tp + fp + fn))


def jaccard(test=None, reference=None, confusion_matrix=None, nan_for_nonexisting=True, **kwargs):
    """TP / (TP + FP + FN)"""

    if confusion_matrix is None:
        confusion_matrix = ConfusionMatrix(test, reference)

    tp, fp, tn, fn = confusion_matrix.get_matrix()
    test_empty, test_full, reference_empty, reference_full = confusion_matrix.get_existence()

    if test_empty and reference_empty:
        if nan_for_nonexisting:
            return float("NaN")
        else:
            return 0.

    return float(tp / (tp + fp + fn))


def precision(test=None, reference=None, confusion_matrix=None, nan_for_nonexisting=True, **kwargs):
    """TP / (TP + FP)"""

    if confusion_matrix is None:
        confusion_matrix = ConfusionMatrix(test, reference)

    tp, fp, tn, fn = confusion_matrix.get_matrix()
    test_empty, test_full, reference_empty, reference_full = confusion_matrix.get_existence()

    if test_empty:
        if nan_for_nonexisting:
            return float("NaN")
        else:
            return 0.

    return float(tp / (tp + fp))


def sensitivity(test=None, reference=None, confusion_matrix=None, nan_for_nonexisting=True, **kwargs):
    """TP / (TP + FN)"""

    if confusion_matrix is None:
        confusion_matrix = ConfusionMatrix(test, reference)

    tp, fp, tn, fn = confusion_matrix.get_matrix()
    test_empty, test_full, reference_empty, reference_full = confusion_matrix.get_existence()

    if reference_empty:
        if nan_for_nonexisting:
            return float("NaN")
        else:
            return 0.

    return float(tp / (tp + fn))


def recall(test=None, reference=None, confusion_matrix=None, nan_for_nonexisting=True, **kwargs):
    """TP / (TP + FN)"""

    return sensitivity(test, reference, confusion_matrix, nan_for_nonexisting, **kwargs)


def specificity(test=None, reference=None, confusion_matrix=None, nan_for_nonexisting=True, **kwargs):
    """TN / (TN + FP)"""

    if confusion_matrix is None:
        confusion_matrix = ConfusionMatrix(test, reference)

    tp, fp, tn, fn = confusion_matrix.get_matrix()
    test_empty, test_full, reference_empty, reference_full = confusion_matrix.get_existence()

    if reference_full:
        if nan_for_nonexisting:
            return float("NaN")
        else:
            return 0.

    return float(tn / (tn + fp))


def accuracy(test=None, reference=None, confusion_matrix=None, **kwargs):
    """(TP + TN) / (TP + FP + FN + TN)"""

    if confusion_matrix is None:
        confusion_matrix = ConfusionMatrix(test, reference)

    tp, fp, tn, fn = confusion_matrix.get_matrix()

    return float((tp + tn) / (tp + fp + tn + fn))


def fscore(test=None, reference=None, confusion_matrix=None, nan_for_nonexisting=True, beta=1., **kwargs):
    """(1 + b^2) * TP / ((1 + b^2) * TP + b^2 * FN + FP)"""

    precision_ = precision(test, reference, confusion_matrix, nan_for_nonexisting)
    recall_ = recall(test, reference, confusion_matrix, nan_for_nonexisting)

    return (1 + beta*beta) * precision_ * recall_ /\
        ((beta*beta * precision_) + recall_)


def false_positive_rate(test=None, reference=None, confusion_matrix=None, nan_for_nonexisting=True, **kwargs):
    """FP / (FP + TN)"""

    return 1 - specificity(test, reference, confusion_matrix, nan_for_nonexisting)


def false_omission_rate(test=None, reference=None, confusion_matrix=None, nan_for_nonexisting=True, **kwargs):
    """FN / (TN + FN)"""

    if confusion_matrix is None:
        confusion_matrix = ConfusionMatrix(test, reference)

    tp, fp, tn, fn = confusion_matrix.get_matrix()
    test_empty, test_full, reference_empty, reference_full = confusion_matrix.get_existence()

    if test_full:
        if nan_for_nonexisting:
            return float("NaN")
        else:
            return 0.

    return float(fn / (fn + tn))


def false_negative_rate(test=None, reference=None, confusion_matrix=None, nan_for_nonexisting=True, **kwargs):
    """FN / (TP + FN)"""

    return 1 - sensitivity(test, reference, confusion_matrix, nan_for_nonexisting)


def true_negative_rate(test=None, reference=None, confusion_matrix=None, nan_for_nonexisting=True, **kwargs):
    """TN / (TN + FP)"""

    return specificity(test, reference, confusion_matrix, nan_for_nonexisting)


def false_discovery_rate(test=None, reference=None, confusion_matrix=None, nan_for_nonexisting=True, **kwargs):
    """FP / (TP + FP)"""

    return 1 - precision(test, reference, confusion_matrix, nan_for_nonexisting)


def negative_predictive_value(test=None, reference=None, confusion_matrix=None, nan_for_nonexisting=True, **kwargs):
    """TN / (TN + FN)"""

    return 1 - false_omission_rate(test, reference, confusion_matrix, nan_for_nonexisting)


def total_positives_test(test=None, reference=None, confusion_matrix=None, **kwargs):
    """TP + FP"""

    if confusion_matrix is None:
        confusion_matrix = ConfusionMatrix(test, reference)

    tp, fp, tn, fn = confusion_matrix.get_matrix()

    return tp + fp


def total_negatives_test(test=None, reference=None, confusion_matrix=None, **kwargs):
    """TN + FN"""

    if confusion_matrix is None:
        confusion_matrix = ConfusionMatrix(test, reference)

    tp, fp, tn, fn = confusion_matrix.get_matrix()

    return tn + fn


def total_positives_reference(test=None, reference=None, confusion_matrix=None, **kwargs):
    """TP + FN"""

    if confusion_matrix is None:
        confusion_matrix = ConfusionMatrix(test, reference)

    tp, fp, tn, fn = confusion_matrix.get_matrix()

    return tp + fn


def total_negatives_reference(test=None, reference=None, confusion_matrix=None, **kwargs):
    """TN + FP"""

    if confusion_matrix is None:
        confusion_matrix = ConfusionMatrix(test, reference)

    tp, fp, tn, fn = confusion_matrix.get_matrix()

    return tn + fp


## ---------- added by Camila
def true_positives(test=None, reference=None, confusion_matrix=None, **kwargs):
    """TP"""

    if confusion_matrix is None:
        confusion_matrix = ConfusionMatrix(test, reference)

    tp, fp, tn, fn = confusion_matrix.get_matrix()

    return tp

def false_positives(test=None, reference=None, confusion_matrix=None, **kwargs):
    """TP"""

    if confusion_matrix is None:
        confusion_matrix = ConfusionMatrix(test, reference)

    tp, fp, tn, fn = confusion_matrix.get_matrix()

    return fp

def true_negatives(test=None, reference=None, confusion_matrix=None, **kwargs):
    """TP"""

    if confusion_matrix is None:
        confusion_matrix = ConfusionMatrix(test, reference)

    tp, fp, tn, fn = confusion_matrix.get_matrix()

    return tn

def false_negatives(test=None, reference=None, confusion_matrix=None, **kwargs):
    """TP"""

    if confusion_matrix is None:
        confusion_matrix = ConfusionMatrix(test, reference)

    tp, fp, tn, fn = confusion_matrix.get_matrix()

    return fn

def relative_volume_error(test=None, reference=None, confusion_matrix=None, **kwargs):
    """abs (P Test - TP Reference) / P Reference"""

    if confusion_matrix is None:
        confusion_matrix = ConfusionMatrix(test, reference)

    tp, fp, tn, fn = confusion_matrix.get_matrix()


    return (abs((tp+fn)-(tp))/(tp+fn))

ALL_METRICS = {
    "False Positive Rate": false_positive_rate,
    "Dice": dice,
    "Jaccard": jaccard,
    # "Hausdorff Distance": hausdorff_distance,
    # "Hausdorff Distance 95": hausdorff_distance_95,
    "Precision": precision,
    "Recall": recall,
    # "Avg. Symmetric Surface Distance": avg_surface_distance_symmetric,
    # "Avg. Surface Distance": avg_surface_distance,
    "Accuracy": accuracy,
    "False Omission Rate": false_omission_rate,
    "Negative Predictive Value": negative_predictive_value,
    "False Negative Rate": false_negative_rate,
    "True Negative Rate": true_negative_rate,
    "False Discovery Rate": false_discovery_rate,
    "Total Positives Test": total_positives_test,
    "Total Negatives Test": total_negatives_test,
    "Total Positives Reference": total_positives_reference,
    "total Negatives Reference": total_negatives_reference,
    ###----------------added by Camila
    "True Positives": true_positives,
    "False Positives": false_positives,
    "True Negatives": true_negatives,
    "False Negatives": false_negatives,
    "Volume Relative Error": relative_volume_error#,
    # "clDice": clDice,
    # "clRecall": clrecall,
    # "clPrecision": clprecision,
    # "clDiceski": skiClDice,
    # "clRecallski": skiClRecall,
    # "clPrecisionski": skiClPrecision
}

# -------------------------------------- evaluator

class Evaluator:
    """Object that holds test and reference segmentations with label information
    and computes a number of metrics on the two. 'labels' must either be an
    iterable of numeric values (or tuples thereof) or a dictionary with string
    names and numeric values.
    """

    default_metrics = [
        "False Positive Rate",
        "Dice",
        "Jaccard",
        "Precision",
        "Recall",
        "Accuracy",
        "False Omission Rate",
        "Negative Predictive Value",
        "False Negative Rate",
        "True Negative Rate",
        "False Discovery Rate",
        "Total Positives Test",
        "Total Positives Reference",
        "True Positives",
        "False Positives",
        "True Negatives",
        "False Negatives",
        "Volume Relative Error"#,
        # "clDice",
        # "clRecall",
        # "clPrecision",
        # "clDiceski",
        # "clRecallski",
        # "clPrecisionski"
    ]

    default_advanced_metrics = [
        #"Hausdorff Distance",
        #"Hausdorff Distance 95",
        #"Avg. Surface Distance",
        #"Avg. Symmetric Surface Distance"
    ]

    def __init__(self,
                 test=None,
                 reference=None,
                 labels=None,
                 metrics=None,
                 advanced_metrics=None,
                 nan_for_nonexisting=True):

        self.test = None
        self.reference = None
        #------- added by Camila
        # self.testcl = None
        # self.referencecl = None
        # -------- end added by Camila
        self.confusion_matrix = ConfusionMatrix()
        self.labels = None
        self.nan_for_nonexisting = nan_for_nonexisting
        self.result = None

        self.metrics = []
        if metrics is None:
            for m in self.default_metrics:
                self.metrics.append(m)
        else:
            for m in metrics:
                self.metrics.append(m)

        self.advanced_metrics = []
        if advanced_metrics is None:
            for m in self.default_advanced_metrics:
                self.advanced_metrics.append(m)
        else:
            for m in advanced_metrics:
                self.advanced_metrics.append(m)

        self.set_reference(reference)
        self.set_test(test)
        if labels is not None:
            self.set_labels(labels)
        else:
            if test is not None and reference is not None:
                self.construct_labels()

    def set_test(self, test):
        """Set the test segmentation."""

        self.test = test

    def set_reference(self, reference):
        """Set the reference segmentation."""

        self.reference = reference

    # ---------- added by Camila

    # def set_testcl(self, testcl):
    #     """Set the test segmentation."""

    #     self.testcl = testcl

    # def set_referencecl(self, referencecl):
    #     """Set the reference segmentation."""

    #     self.referencecl = referencecl

    # ---------- end added by Camila

    def set_labels(self, labels):
        """Set the labels.
        :param labels= may be a dictionary (int->str), a set (of ints), a tuple (of ints) or a list (of ints). Labels
        will only have names if you pass a dictionary"""

        if isinstance(labels, dict):
            self.labels = collections.OrderedDict(labels)
        elif isinstance(labels, set):
            self.labels = list(labels)
        elif isinstance(labels, np.ndarray):
            self.labels = [i for i in labels]
        elif isinstance(labels, (list, tuple)):
            self.labels = labels
        else:
            raise TypeError("Can only handle dict, list, tuple, set & numpy array, but input is of type {}".format(type(labels)))

    def construct_labels(self):
        """Construct label set from unique entries in segmentations."""

        if self.test is None and self.reference is None:
            raise ValueError("No test or reference segmentations.")
        elif self.test is None:
            labels = np.unique(self.reference)
        else:
            labels = np.union1d(np.unique(self.test),
                                np.unique(self.reference))
        self.labels = list(map(lambda x: int(x), labels))

    def set_metrics(self, metrics):
        """Set evaluation metrics"""

        if isinstance(metrics, set):
            self.metrics = list(metrics)
        elif isinstance(metrics, (list, tuple, np.ndarray)):
            self.metrics = metrics
        else:
            raise TypeError("Can only handle list, tuple, set & numpy array, but input is of type {}".format(type(metrics)))

    def add_metric(self, metric):

        if metric not in self.metrics:
            self.metrics.append(metric)

    def evaluate(self, test=None, reference=None, advanced=False, **metric_kwargs):
        """Compute metrics for segmentations."""
        if test is not None:
            self.set_test(test)

        if reference is not None:
            self.set_reference(reference)

        if self.test is None or self.reference is None:
            raise ValueError("Need both test and reference segmentations.")

        if self.labels is None:
            self.construct_labels()

        self.metrics.sort()

        # get functions for evaluation
        # somewhat convoluted, but allows users to define additonal metrics
        # on the fly, e.g. inside an IPython console
        _funcs = {m: ALL_METRICS[m] for m in self.metrics + self.advanced_metrics}
        frames = inspect.getouterframes(inspect.currentframe())
        for metric in self.metrics:
            for f in frames:
                if metric in f[0].f_locals:
                    _funcs[metric] = f[0].f_locals[metric]
                    break
            else:
                if metric in _funcs:
                    continue
                else:
                    raise NotImplementedError(
                        "Metric {} not implemented.".format(metric))

        # get results
        self.result = OrderedDict()

        eval_metrics = self.metrics
        if advanced:
            eval_metrics += self.advanced_metrics

        if isinstance(self.labels, dict):

            for label, name in self.labels.items():
                k = str(name)
                self.result[k] = OrderedDict()
                if not hasattr(label, "__iter__"):
                    self.confusion_matrix.set_test(self.test == label)
                    self.confusion_matrix.set_reference(self.reference == label)
                    # --- added by Camila
                    # self.confusion_matrix.set_testcl(self.testcl == label)
                    # self.confusion_matrix.set_referencecl(self.referencecl == label)
                    # self.confusion_matrix.set_skirefcl(skeletonize(self.reference)==label)
                    # self.confusion_matrix.set_skitestcl(skeletonize(self.test)==label)
                    ## ------------------
                else:
                    current_test = 0
                    current_reference = 0
                    for l in label:
                        current_test += (self.test == l)
                        current_reference += (self.reference == l)
                    self.confusion_matrix.set_test(current_test)
                    self.confusion_matrix.set_reference(current_reference)
                for metric in eval_metrics:
                    self.result[k][metric] = _funcs[metric](confusion_matrix=self.confusion_matrix,
                                                               nan_for_nonexisting=self.nan_for_nonexisting,
                                                               **metric_kwargs)

        else:

            for i, l in enumerate(self.labels):
                k = str(l)
                self.result[k] = OrderedDict()
                self.confusion_matrix.set_test(self.test == l)
                self.confusion_matrix.set_reference(self.reference == l)
                # --- added by Camila
                # self.confusion_matrix.set_testcl(self.testcl == l)
                # self.confusion_matrix.set_referencecl(self.referencecl == l)
                # self.confusion_matrix.set_skirefcl(skeletonize(self.reference)==l)
                # self.confusion_matrix.set_skitestcl(skeletonize(self.test)==l)
                ## ------------------
                for metric in eval_metrics:
                    self.result[k][metric] = _funcs[metric](confusion_matrix=self.confusion_matrix,
                                                            nan_for_nonexisting=self.nan_for_nonexisting,
                                                            **metric_kwargs)

        return self.result

    def to_dict(self):

        if self.result is None:
            self.evaluate()
        return self.result

    def to_array(self):
        """Return result as numpy array (labels x metrics)."""

        if self.result is None:
            self.evaluate

        result_metrics = sorted(self.result[list(self.result.keys())[0]].keys())

        a = np.zeros((len(self.labels), len(result_metrics)), dtype=np.float32)

        if isinstance(self.labels, dict):
            for i, label in enumerate(self.labels.keys()):
                for j, metric in enumerate(result_metrics):
                    a[i][j] = self.result[self.labels[label]][metric]
        else:
            for i, label in enumerate(self.labels):
                for j, metric in enumerate(result_metrics):
                    a[i][j] = self.result[label][metric]

        return a


class NiftiEvaluator(Evaluator):

    def __init__(self, *args, **kwargs):

        self.test_nifti = None
        self.reference_nifti = None
        self.test_nifti_cl = None
        self.reference_nifti_cl = None
        super(NiftiEvaluator, self).__init__(*args, **kwargs)

    def set_test(self, test):
        """Set the test segmentation."""

        if test is not None:
            self.test_nifti = sitk.ReadImage(test)
            super(NiftiEvaluator, self).set_test(sitk.GetArrayFromImage(self.test_nifti))
        else:
            self.test_nifti = None
            super(NiftiEvaluator, self).set_test(test)

    def set_reference(self, reference):
        """Set the reference segmentation."""

        if reference is not None:
            self.reference_nifti = sitk.ReadImage(reference)
            super(NiftiEvaluator, self).set_reference(sitk.GetArrayFromImage(self.reference_nifti))
        else:
            self.reference_nifti = None
            super(NiftiEvaluator, self).set_reference(reference)

    # -------- added by Camila

    # def computeThinningSITK(self, inputImage):
    #     thinningFilter = sitk.BinaryThinningImageFilter()
    #     thinnedImage = thinningFilter.Execute(inputImage)
    #     return thinnedImage

    # def set_testcl(self):
    #     """Set the test segmentation."""

    #     if self.test_nifti is not None:            
    #         self.test_nifticl = self.computeThinningSITK(self.test_nifti)
    #         super(NiftiEvaluator, self).set_testcl(sitk.GetArrayFromImage(self.test_nifticl))
    #     else:
    #         self.test_nifticl = None
    #         super(NiftiEvaluator, self).set_testcl(None)

    # def set_referencecl(self):
    #     """Set the reference segmentation."""

    #     if self.reference_nifti is not None:
    #         self.reference_nifticl = self.computeThinningSITK(self.reference_nifti)
    #         super(NiftiEvaluator, self).set_referencecl(sitk.GetArrayFromImage(self.reference_nifticl))
    #     else:
    #         self.reference_nifticl = None
    #         super(NiftiEvaluator, self).set_referencecl(None)

    # ----------- end added by Camila

    def evaluate(self, test=None, reference=None, voxel_spacing=None, **metric_kwargs):

        if voxel_spacing is None:
            voxel_spacing = np.array(self.test_nifti.GetSpacing())[::-1]
            metric_kwargs["voxel_spacing"] = voxel_spacing

        return super(NiftiEvaluator, self).evaluate(test, reference, **metric_kwargs)


def run_evaluation(args):
    test, ref, evaluator, metric_kwargs = args
    # evaluate
    evaluator.set_test(test)
    evaluator.set_reference(ref)
    #------ added by Camila
    # evaluator.set_testcl()
    # evaluator.set_referencecl()
    #---------- end added by Camila
    if evaluator.labels is None:
        evaluator.construct_labels()
    current_scores = evaluator.evaluate(**metric_kwargs)
    if type(test) == str:
        current_scores["test"] = test
    if type(ref) == str:
        current_scores["reference"] = ref
    return current_scores


def aggregate_scores(test_ref_pairs,
                     evaluator=NiftiEvaluator,
                     labels=None,
                     nanmean=True,
                     json_output_file=None,
                     json_name="",
                     json_description="",
                     json_author="Camila",
                     json_task="",
                     num_threads=2,
                     **metric_kwargs):
    """
    test = predicted image
    :param test_ref_pairs:
    :param evaluator:
    :param labels: must be a dict of int-> str or a list of int
    :param nanmean:
    :param json_output_file:
    :param json_name:
    :param json_description:
    :param json_author:
    :param json_task:
    :param metric_kwargs:
    :return:
    """

    if type(evaluator) == type:
        evaluator = evaluator()

    if labels is not None:
        evaluator.set_labels(labels)

    all_scores = OrderedDict()
    all_scores["all"] = []
    all_scores["mean"] = OrderedDict()
    all_scores["sum"] = OrderedDict()

    test = [i[0] for i in test_ref_pairs]
    ref = [i[1] for i in test_ref_pairs]
    p = Pool(num_threads)
    all_res = p.map(run_evaluation, zip(test, ref, [evaluator]*len(ref), [metric_kwargs]*len(ref)))
    p.close()
    p.join()

    for i in range(len(all_res)):
        all_scores["all"].append(all_res[i])

        # append score list for mean
        for label, score_dict in all_res[i].items():
            if label in ("test", "reference"):
                continue
            if label not in all_scores["mean"]:
                all_scores["mean"][label] = OrderedDict()
            # ------- added by Camila
            if label not in all_scores["sum"]:
                all_scores["sum"][label] = OrderedDict()
            # -----------
            for score, value in score_dict.items():
                if score not in all_scores["mean"][label]:
                    all_scores["mean"][label][score] = []
                # ------- added by Camila
                if score not in all_scores["sum"][label] and score in ["True Positives", "False Positives", "True Negatives", "False Negatives", "clDice"]:
                    all_scores["sum"][label][score] = []
                # ----------    
                all_scores["mean"][label][score].append(value)
                # ------- added by Camila
                if score in ["True Positives", "False Positives", "True Negatives", "False Negatives", "clDice"]:
                    all_scores["sum"][label][score].append(value)
                # ---------------

    for label in all_scores["mean"]:
        for score in all_scores["mean"][label]:
            if nanmean:
                all_scores["mean"][label][score] = float(np.nanmean(all_scores["mean"][label][score]))
            else:
                all_scores["mean"][label][score] = float(np.mean(all_scores["mean"][label][score]))
    
    # ------- added by Camila
    for label in all_scores["sum"]:
        for score in all_scores["sum"][label]:
            all_scores["sum"][label][score] = float(np.sum(all_scores["sum"][label][score]))
    #for label in all_scores["sum"]:
        all_scores["sum"][label]["Sensitivity"] = float(all_scores["sum"][label]["True Positives"] / (all_scores["sum"][label]["True Positives"] + all_scores["sum"][label]["False Negatives"]))
        all_scores["sum"][label]["Specificity"] = float(all_scores["sum"][label]["True Negatives"] / (all_scores["sum"][label]["True Negatives"] + all_scores["sum"][label]["False Positives"]))
        all_scores["sum"][label]["Precision"] = float(all_scores["sum"][label]["True Positives"] / (all_scores["sum"][label]["True Positives"] + all_scores["sum"][label]["False Positives"]))
        all_scores["sum"][label]["Accuracy"] = float((all_scores["sum"][label]["True Positives"] + all_scores["sum"][label]["True Negatives"])/ (all_scores["sum"][label]["True Positives"] + all_scores["sum"][label]["False Negatives"] + all_scores["sum"][label]["False Positives"] + all_scores["sum"][label]["True Negatives"]))
        all_scores["sum"][label]["DICE"] = float(2.*(all_scores["sum"][label]["True Positives"])/ (2*all_scores["sum"][label]["True Positives"] + all_scores["sum"][label]["False Negatives"] + all_scores["sum"][label]["False Positives"]))
        #all_scores["sum"][label]["clDice"] = float(all_scores["sum"][label]["clDice"]/len(all_res))
    # ------------------

    # save to file if desired
    # we create a hopefully unique id by hashing the entire output dictionary
    if json_output_file is not None:
        json_dict = OrderedDict()
        json_dict["name"] = json_name
        json_dict["description"] = json_description
        timestamp = datetime.today()
        json_dict["timestamp"] = str(timestamp)
        json_dict["task"] = json_task
        json_dict["author"] = json_author
        json_dict["results"] = all_scores
        json_dict["id"] = hashlib.md5(json.dumps(json_dict).encode("utf-8")).hexdigest()[:12]
        writejson(json_dict, json_output_file)


    return all_scores


def aggregate_scores_for_experiment(score_file,
                                    labels=None,
                                    metrics=Evaluator.default_metrics,
                                    nanmean=True,
                                    json_output_file=None,
                                    json_name="",
                                    json_description="",
                                    json_author="Fabian",
                                    json_task=""):

    scores = np.load(score_file)
    scores_mean = scores.mean(0)
    if labels is None:
        labels = list(map(str, range(scores.shape[1])))

    results = []
    results_mean = OrderedDict()
    for i in range(scores.shape[0]):
        results.append(OrderedDict())
        for l, label in enumerate(labels):
            results[-1][label] = OrderedDict()
            results_mean[label] = OrderedDict()
            for m, metric in enumerate(metrics):
                results[-1][label][metric] = float(scores[i][l][m])
                results_mean[label][metric] = float(scores_mean[l][m])

    json_dict = OrderedDict()
    json_dict["name"] = json_name
    json_dict["description"] = json_description
    timestamp = datetime.today()
    json_dict["timestamp"] = str(timestamp)
    json_dict["task"] = json_task
    json_dict["author"] = json_author
    json_dict["results"] = {"all": results, "mean": results_mean}
    json_dict["id"] = hashlib.md5(json.dumps(json_dict).encode("utf-8")).hexdigest()[:12]
    if json_output_file is not None:
        json_output_file = open(json_output_file, "w")
        json.dump(json_dict, json_output_file, indent=4, separators=(",", ": "))
        json_output_file.close()

    return json_dict


def evaluate_folder(folder_with_gts: str, folder_with_predictions: str, labels: tuple, **metric_kwargs):
    """
    writes a summary.json to folder_with_predictions
    :param folder_with_gts: folder where the ground truth segmentations are saved. Must be nifti files.
    :param folder_with_predictions: folder where the predicted segmentations are saved. Must be nifti files.
    :param labels: tuple of int with the labels in the dataset. For example (0, 1, 2, 3) for Task001_BrainTumour.
    :return:
    """
    l = lambda x, y: y
    #res = [l(folder_with_gts, i) for i in os.listdir(folder_with_gts) if os.path.isfile(os.path.join(folder_with_gts, i))
    #        and (i.endswith(".nii.gz"))].sort()
    files_gt = [l(folder_with_gts, i) for i in os.listdir(folder_with_gts) if os.path.isfile(os.path.join(folder_with_gts, i))
            and (i.endswith(".nii.gz"))]
    files_gt.sort()
    files_pred = [l(folder_with_predictions, i) for i in os.listdir(folder_with_predictions) if os.path.isfile(os.path.join(folder_with_predictions, i))
            and (i.endswith(".nii.gz"))]
    files_pred.sort()
    assert all([i in files_pred for i in files_gt]), "files missing in folder_with_predictions"
    assert all([i in files_gt for i in files_pred]), "files missing in folder_with_gts"
    test_ref_pairs = [(os.path.join(folder_with_predictions, i), os.path.join(folder_with_gts, i)) for i in files_pred]
    res = aggregate_scores(test_ref_pairs, json_output_file=os.path.join(folder_with_predictions, "summary.json"),
                           num_threads=1, labels=labels, **metric_kwargs)
    return res


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-ref', required=True, type=str, help="Folder containing the reference segmentations in nifti "
                                                              "format.")
    parser.add_argument('-pred', required=True, type=str, help="Folder containing the predicted segmentations in nifti "
                                                               "format. File names must match between the folders!")
    parser.add_argument('-l', nargs='+', type=int, required=True, help="List of labels to evaluate (-l 1)")
    args = parser.parse_args()
    return evaluate_folder(args.ref, args.pred, args.l)


if __name__ == "__main__":
    main()
