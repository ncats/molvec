package tripod.molvec.algo;

import org.junit.Ignore;
import org.junit.experimental.ParallelComputer;
import org.junit.runner.Computer;
import org.junit.runner.JUnitCore;
import org.junit.runner.Result;
import org.junit.runner.Runner;
import org.junit.runners.ParentRunner;
import org.junit.runners.model.InitializationError;
import org.junit.runners.model.RunnerBuilder;
import org.junit.runners.model.RunnerScheduler;

import java.time.Duration;
import java.time.Period;
import java.util.Objects;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Created by katzelda on 3/12/19.
 */
@Ignore
public class ParallelRegressionTest2 {

    public static void main(String[] args){
        Class[] cls = new Class[]{RegressionTest2.class};


        Result result = JUnitCore.runClasses(MyParallelComputer.methods(Executors.newFixedThreadPool(3)), cls);
        result.getFailures().forEach(System.err::println);
        System.out.printf("Tests Complete: %d total %d failed %d Ignored%nRun Time = %s%n",
                result.getRunCount(), result.getFailureCount(), result.getIgnoreCount(),
                Duration.ofMillis(result.getRunTime()));
    }

    public static class MyParallelComputer extends Computer {
        private final boolean classes;

        private final boolean methods;

        private final ExecutorService executorService;

        public MyParallelComputer(ExecutorService executorService, boolean classes, boolean methods) {
            this.classes = classes;
            this.methods = methods;
            this.executorService = Objects.requireNonNull(executorService);
        }

        public static Computer classes(ExecutorService executorService) {
            return new MyParallelComputer(executorService, true, false);
        }

        public static Computer methods(ExecutorService executorService) {
            return new MyParallelComputer(executorService, false, true);
        }

        private Runner parallelize(Runner runner) {
            if (runner instanceof ParentRunner) {
                ((ParentRunner<?>) runner).setScheduler(new RunnerScheduler() {

                    public void schedule(Runnable childStatement) {
                        executorService.submit(childStatement);
                    }

                    public void finished() {
                        try {
                            executorService.shutdown();
                            executorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
                        } catch (InterruptedException e) {
                            e.printStackTrace(System.err);
                        }
                    }
                });
            }
            return runner;
        }

        @Override
        public Runner getSuite(RunnerBuilder builder, java.lang.Class<?>[] classes)
                throws InitializationError {
            Runner suite = super.getSuite(builder, classes);
            return this.classes ? parallelize(suite) : suite;
        }

        @Override
        protected Runner getRunner(RunnerBuilder builder, Class<?> testClass)
                throws Throwable {
            Runner runner = super.getRunner(builder, testClass);
            return methods ? parallelize(runner) : runner;
        }
    }
}
