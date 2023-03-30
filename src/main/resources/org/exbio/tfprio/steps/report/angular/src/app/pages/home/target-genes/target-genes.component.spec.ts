import { ComponentFixture, TestBed } from '@angular/core/testing';

import { TargetGenesComponent } from './target-genes.component';

describe('TargetGenesComponent', () => {
  let component: TargetGenesComponent;
  let fixture: ComponentFixture<TargetGenesComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ TargetGenesComponent ]
    })
    .compileComponents();

    fixture = TestBed.createComponent(TargetGenesComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
