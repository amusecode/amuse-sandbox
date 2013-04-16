import nl.esciencecenter.amuse.distributed.DistributedAmuse;

public class Code implements CodeInterface {

    public Code() {
        new DistributedAmuse();

    }
    

    @Override
    public int initialize_code() {
        // TODO Auto-generated method stub
        return 0;
    }
    
    @Override
    public int get_port() {
        return 56453;
    }

    @Override
    public int new_resource(int[] index_of_the_resource, String[] name, String[] hostname, String[] username,
            String[] scheduler, String[] amuse_dir, int count) {
        // TODO Auto-generated method stub
        return 0;
    }
    
    @Override
    public int delete_resource(int[] index_of_the_resource, int npoints) {
        // TODO Auto-generated method stub
        return 0;
    }
    
    @Override
    public int new_reservation(int[] reservation_id, String[] resource_name, int[] node_count, int[] time,
            String[] node_label) {
        // TODO Auto-generated method stub
        return 0;
    }
    
    @Override
    public int delete_reservation(int[] reservation_id, int count) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int submit_pickled_function_job(int[] job_id, String[] function, String[] arguments, String[] node_label,
            int count) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int submit_script_job(int[] job_id, String[] script, String[] arguments, String[] code_dir,
            String[] node_label, int[] re_use_code_files, int count) {
        // TODO Auto-generated method stub
        return 0;
    }
    
    @Override
    public int get_pickled_function_job_result(int[] job_id, String[] result, int count) {
        // TODO Auto-generated method stub
        for (int i = 0; i < count; i++) {
            result[i] = "lots of pickles!";
        }
        
        return 0;
    }

    @Override
    public int wait_for_jobs() {
        // TODO Auto-generated method stub
        return 0;
    }
   
    @Override
    public int recommit_parameters() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int wait_for_reservations() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int commit_parameters() {
        // TODO Auto-generated method stub
        return 0;
    }
}
