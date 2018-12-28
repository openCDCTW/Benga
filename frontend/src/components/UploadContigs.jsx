import React from 'react';
import ReactDOM from 'react-dom';
import Navigation from './Navigation.jsx';
import DropzoneComponent from 'react-dropzone-component';
import { Link } from 'react-router-dom';
import Options from './Options.jsx';
import Button from '@material-ui/core/Button';
import { withStyles } from '@material-ui/core/styles';
import CloudUploadIcon from '@material-ui/icons/CloudUpload';
import Icon from '@material-ui/core/Icon';
import DeleteIcon from '@material-ui/icons/Delete';
import ReplyIcon from '@material-ui/icons/Reply';
import DownloadIcon from '@material-ui/icons/CloudDownload';
import SearchBar from './SearchBar.jsx';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import green from '@material-ui/core/colors/green';
// import { Scrollbars } from 'react-custom-scrollbars';

const styles = theme => ({
    cssRoot:{
        color: theme.palette.getContrastText(green[600]),
        backgroundColor: green[500],
        '&:hover': {
            backgroundColor:green[600],
        },
    }
})

class Upload_contigs extends React.Component {

    constructor(props) {

        super(props);

        fetch('api/profiling/upload/', {method:'POST'})
        .then(function(res){
           return res.json();
        }).then(function(batch){
           return getID(batch);
        });

        var getID=function(data){
            window.batchid = data.id;
        };


        // TODO: poor performnce issue with everytime load the component will
        // fetch API once.

        window.databaseName = "";
        window.fileName = [];

        this.state = {
            to: "/",
            upload_confirm: false,
            switch: false,
        };

        this.djsConfig = {
            dictDefaultMessage:"Drop files or click to upload contigs",
            addRemoveLinks: true,
            acceptedFiles: ".fasta,.fa,.fna",
            autoProcessQueue: false,
            parallelUploads: 200,
            init:function(){
                this.on("sending", function(file, xhr, formData){
                    formData.append("batch_id", window.batchid);
                });
                this.on("success", function(file){
                    window.fileName.push(file.name);
                    file._removeLink.remove();
                    delete file._removeLink;
                });
            }
        }

        this.componentConfig = {
            iconFiletypes: ['.fasta','.fna','.fa'],
            showFiletypeIcon: true,
            postUrl: 'api/profiling/sequence/'
        };

        this.dropzone = null;
    }

    handlePost() {
        
        let fileCheck = this.dropzone.files.length;

        if(fileCheck < 5){
            alert('Please upload at least 5 files');
            return ;
        }

        if(window.databaseName == ""){
            alert('Please choose a database !');
            return ;
        }


        this.dropzone.processQueue();
        this.setState(state => ({ upload_confirm: true , to: '/profile_view' ,
            switch: true }));
    }

    submit(){


        if (this.state.upload_confirm == false){
            alert('Please upload files first! (At least 5 files)');
            return ;
        }

        var scheme = {};
        scheme.occurrence = "95";
        scheme.database = window.databaseName;
        scheme.id = window.batchid;
        fetch('api/profiling/profiling-tree/', {
            method:'POST',
            headers: new Headers({'content-type': 'application/json'}),
            body: JSON.stringify(scheme)
        });

    }

    remove(){
        
        this.dropzone.removeAllFiles();
        this.setState(state => ({ to: '/', upload_confirm: false, switch: false }));
        fetch('api/profiling/upload/', {method:'POST'})
        .then(function(res){
           return res.json();
        }).then(function(batch){
           return getID(batch);
        });

        var getID=function(data){
            window.batchid = data.id;
        };
    }

    query(){

        fetch('api/profiling/profile/' + this.state.querybyID, { method:'GET'})
        .then(response => response.json())
        .catch(error => alert("Data not found"));

        fetch('api/profiling/profile/' + this.state.querybyID, { method:'GET'})
            .then(response => response.json())
            .then(result => this.setState(state => ({
                profile_result_all: result.file,
                profile_result_zip: result.zip,})));

        fetch('api/dendrogram/dendrogram/' + this.state.querybyID, { method: 'GET'})
            .then(response => response.json())
            .then(result => this.setState(state => ({
                png_file: result.png_file, 
                pdf_file: result.pdf_file,
                svg_file: result.svg_file, 
                emf_file: result.emf_file, 
                newick_file: result.newick_file, })));
    }

    upload_samaple(){
        
    }

    back(){
        window.location.reload();
    }

    
    render() {

        const config = this.componentConfig;
        const djsConfig = this.djsConfig;
        const eventHandlers = {
            init: dz => this.dropzone = dz,}
        const { classes } = this.props;

        if(this.state.profile_result_all == undefined){
            return (
            <div>
                <Navigation value={0}/>
                <br />
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                    <SearchBar
                        onChange = {(value) => this.setState({ querybyID: value })}
                        onRequestSearch={this.query.bind(this)}
                        placeholder = {"Input batch ID here to query data"}
                        style = {{
                            width: '90%',
                            margin: '0 auto',
                        }}
                    />
                </div>
                <br />
                <div style={{ width:'97%', display:'flex', justifyContent:'flex-end', 
                alignItems:'flex-end'}}>
                    <Button variant="contained" color="secondary" onClick={this.remove.bind(this)}>
                            Remove all files
                            &nbsp;&nbsp;
                            <DeleteIcon />
                    </Button>
                </div>
                <br />
                <div style = {{ display:'flex', justifyContent:'center', alignItems:'center' }}>
                    <DropzoneComponent config={config} eventHandlers={eventHandlers} 
                        djsConfig={djsConfig} />
                </div>
                <br />
                <br />
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                    <font> SUGGESTION: Please use &nbsp;
                        <a href="http://cab.spbu.ru/software/spades/" target="_blank">SPAdes</a>
                        &nbsp; to assembly before upload. 
                    </font>
                </div>
                <br />
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                    <Options switch={this.state.switch} />
                </div>
                <br />
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                    <Button variant="contained" className ={classes.cssRoot} 
                     onClick={this.handlePost.bind(this)}>
                        Upload
                        &nbsp;&nbsp;
                        <CloudUploadIcon />
                    </Button>
                    &nbsp;&nbsp;&nbsp;&nbsp;
                    <Link to={this.state.to} style={{ textDecoration:'none' }}>
                        <Button variant="contained" color="primary" onClick={this.submit.bind(this)}>
                            profiling
                            &nbsp;&nbsp;
                            <Icon>send</Icon>
                        </Button>
                    </Link>
                </div>
                <br />
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                    <Button variant="contained" color="default" 
                     onClick={this.upload_samaple.bind(this)}>
                        Upload sample data
                        &nbsp;&nbsp;
                        <CloudUploadIcon />
                    </Button>
                </div>
                <br />
                <br />
                <br />
                <br />
            </div>
            );
        }else{
            return (
                    <div>
                        <Paper square>
                            <Tabs centered>
                                <Tab label=" "/>
                            </Tabs>
                        </Paper>
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <a download href={this.state.profile_result_all} 
                             style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Download profiles (.tsv)
                                &nbsp;&nbsp;
                                <DownloadIcon />
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.profile_result_zip} 
                             style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Download profiles (.zip)
                                &nbsp;&nbsp;
                                <DownloadIcon />
                                </Button>
                            </a>
                        </div>
                        <br />
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <img src={this.state.svg_file} />
                        </div>
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <font size="4">Download Dendrogram</font>
                        </div>
                        <br />
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            
                            <a download href={this.state.png_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Png 
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.pdf_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Pdf
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.svg_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Svg
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.emf_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                emf
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.newick_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                newick
                                </Button>
                            </a>
                        </div>
                        <br />
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center' }}>
                            <Link to="/" style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default" onClick={this.back.bind(this)}>
                                    <ReplyIcon />
                                    &nbsp;&nbsp;
                                    Back
                                </Button>
                            </Link>
                        </div>
                        <br />
                        <br />
                    </div>
                );
        }
        
    }
}


export default withStyles(styles)(Upload_contigs);