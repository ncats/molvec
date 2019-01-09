package tripod.molvec;


import java.util.Optional;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Supplier;

/**
 * Memoized supplier. Caches the result of the supplier
 * to be used after. Useful for expensive calls.
 * 
 * @author peryeata
 * @param <T>
 */
public class CachedSupplier<T> implements Supplier<T>, Callable<T>{
	private static AtomicLong generatedVersion= new AtomicLong();
	
	
	

	/**
	 * Flag to signal all {{@link ix.core.util.CachedSupplier} instances
	 * to regenerate from their suppliers on the next call.
	 */
	public static void resetAllCaches(){
		CachedSupplier.generatedVersion.incrementAndGet();
	}

	private final Supplier<T> c;
	private T cache;
	private boolean run=false;
	private long generatedWithVersion;

	public CachedSupplier(final Supplier<T> c){
		this.c=c;
	}

	/**
	 * Delegates to {@link #get()}
	 */
	@Override
	public T call() throws Exception{
		return get();
	}
		
	
	@Override
	public T get() {
		if(hasRun()) {
			return this.cache;
		}else{
			synchronized(this){
				if(hasRun()){
					return this.cache;
				}
				this.generatedWithVersion=CachedSupplier.generatedVersion.get();
				this.cache=directCall();
				this.run=true;
				return this.cache;
			}
		}
	}
	
	protected T directCall(){
		
		return this.c.get();
	}


	/**
	 * An explicitly synchronized form of {@link #get()}
	 * @return
	 */
	public synchronized T getSync() {
		return get();
	}
	
	public boolean hasRun(){
		return this.run && this.generatedWithVersion==CachedSupplier.generatedVersion.get();
	}


	/**
	 * Flag to signal this instance to recalculate from its
	 * supplier on next call.
	 */
	public void resetCache(){
		this.run=false;
	}

	public static <T> CachedSupplier<T> of(final Supplier<T> supplier){
		return new CachedSupplier<T>(supplier);
	}

	/**
	 * Wrap the provided callable as a cached supplier
	 * @param callable
	 * @return
	 */
	public static <T> CachedSupplier<T> ofCallable(final Callable<T> callable){
		return of(()->{
			try{
				return callable.call();
			}catch(final Exception e){
				throw new IllegalStateException(e);
			}
		});
	}
	
	public static <T> CachedThrowingSupplier<T> ofThrowing(final Callable<T> callable){
		return new CachedThrowingSupplier<T>(()->{
			try{
				return callable.call();
			}catch(final Exception e){
				throw new IllegalStateException(e);
			}
		});
	}
	
	/**
	 * An extension of a {@link CachedSupplier} which will catch any
	 * throwable thrown during the initial {@link Supplier#get()} call,
	 * and cache it as well, returning <code>null</code> for the value
	 * cache. Calling {@link #getThrown()} will return an {@link Optional}
	 * of a {@link Throwable}, which is empty if there was nothing 
	 * thrown during the execution.
	 * @author peryeata
	 *
	 * @param <T>
	 */
	public static class CachedThrowingSupplier<T> extends CachedSupplier<T>{

		public Throwable thrown=null;
		
		public CachedThrowingSupplier(Supplier<T> c) {
			super(c);
		}
		
		@Override
		protected T directCall(){
			try{
				return super.directCall();
			}catch(Throwable e){
				setThrown(e);
				return null;
			}
		}
		
		private void setThrown(Throwable t){
			this.thrown=t;
		}
		
		/**
		 * Calls the supplier (if necessary), and returns an {@link Optional}
		 * of anything thrown by that supplier.
		 * @return
		 */
		public Optional<Throwable> getThrown(){
			this.get();
			return Optional.ofNullable(thrown);
		}
		
	}
}